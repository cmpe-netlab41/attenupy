"""Study designs: samplers that produce :class:`~attenupy.spec.SimulationSpec`.

A *sampler* is anything iterable that yields specs -- a plain list works too.
Subclass :class:`Sampler` (implement ``__iter__``) to get ``+`` chaining and
``.sample()`` for free. Built in:

* :class:`RandomSampler`  -- N random mixtures per beam, random thickness each
  (the paper's training design).
* :class:`GridSampler`    -- beams x mixtures x thicknesses (full product).
* :class:`ThicknessSweep` -- given mixture(s), every thickness in a list
  (the paper's held-out test curves).
* :class:`CSVSampler`     -- one or more specs per CSV row.

The composition logic is pluggable too:

* :class:`RandomMixtures` draws mixtures: how many constituents, which
  combination, which polymer fraction, how the rest is split.
* A *splitter* (:class:`Dirichlet`, :class:`EqualSplit`, or any object with
  ``split(names, total, rng) -> list[float]``) divides the non-polymer mass.

All randomness comes from a ``random.Random(seed)`` created per iteration, so a
sampler with a fixed ``seed`` yields the same design every time it is iterated.
"""

import csv
import itertools
import logging
import math
import random
from abc import ABC, abstractmethod

from .data import POLYMER_FORMULA
from .spec import Mixture, SimulationSpec

log = logging.getLogger(__name__)


class Sampler(ABC):
    """Base class: iterate to get specs."""

    @abstractmethod
    def __iter__(self):
        """Yield :class:`SimulationSpec` objects."""

    def sample(self):
        return list(self)

    def __add__(self, other):
        return Chain(self, other)


class Chain(Sampler):
    """Concatenation of several samplers (what ``a + b`` returns)."""

    def __init__(self, *parts):
        self.parts = parts

    def __iter__(self):
        for part in self.parts:
            yield from part


# ---------------------------------------------------------------------------
# Splitting the non-polymer mass
# ---------------------------------------------------------------------------
class Dirichlet:
    """Split mass by a Dirichlet draw, optionally with a guaranteed floor.

    ``Dirichlet()`` is the flat Dirichlet of the paper. ``alpha`` may be a
    number (all constituents) or ``{name: alpha}`` (others get
    ``default_alpha``), e.g. ``Dirichlet({"B": 2.5, "Ba": 2.5}, min_share=0.1)``
    biases B/Ba toward larger shares while keeping the draw random.
    ``min_share`` reserves that fraction of the split mass for *each* constituent.
    """

    def __init__(self, alpha=1.0, *, default_alpha=1.0, min_share=0.0):
        self.alpha = alpha
        self.default_alpha = default_alpha
        self.min_share = min_share

    def _alpha(self, name):
        if isinstance(self.alpha, (int, float)):
            return float(self.alpha)
        return float(self.alpha.get(name, self.default_alpha))

    def split(self, names, total, rng):
        k = len(names)
        if self.min_share * k >= 1.0:
            raise ValueError(f"min_share={self.min_share} is too large for {k} constituents")
        draws = [rng.gammavariate(self._alpha(n), 1.0) for n in names]
        s = sum(draws)
        base = self.min_share * total
        free = total - k * base
        return [base + d / s * free for d in draws]


class EqualSplit:
    """Every constituent gets the same share."""

    def split(self, names, total, rng):
        return [total / len(names)] * len(names)


def round_preserving_sum(values, decimals):
    """Round to ``decimals`` so the rounded values still sum to the rounded
    total, and none rounds to zero (each keeps at least one unit)."""
    scale = 10 ** decimals
    total_units = round(sum(values) * scale)
    k = len(values)
    if total_units < k:
        raise ValueError(f"Cannot give {k} constituents a non-zero share of "
                         f"{sum(values)} at {decimals} decimals")
    rest = total_units - k
    raw = [max(v * scale - 1.0, 0.0) for v in values]
    s = sum(raw)
    shares = [r / s * rest for r in raw] if s > 0 else [rest / k] * k
    units = [math.floor(x) for x in shares]
    by_remainder = sorted(range(k), key=lambda i: shares[i] - units[i], reverse=True)
    for i in by_remainder[: rest - sum(units)]:
        units[i] += 1
    return [round((u + 1) / scale, decimals) for u in units]


# ---------------------------------------------------------------------------
# Random mixtures
# ---------------------------------------------------------------------------
def _normalise_weights(n_constituents):
    if isinstance(n_constituents, int):
        return {n_constituents: 1.0}
    if hasattr(n_constituents, "items"):
        weights = {int(k): float(w) for k, w in n_constituents.items() if w > 0}
    else:
        weights = {int(k): 1.0 for k in n_constituents}
    if not weights:
        raise ValueError("n_constituents is empty")
    total = sum(weights.values())
    return {k: w / total for k, w in weights.items()}


def _allocate(count, weights):
    """Largest-remainder split of ``count`` according to ``weights``."""
    exact = {k: count * w for k, w in weights.items()}
    alloc = {k: math.floor(v) for k, v in exact.items()}
    for k in sorted(exact, key=lambda k: exact[k] - alloc[k], reverse=True)[: count - sum(alloc.values())]:
        alloc[k] += 1
    return alloc


class RandomMixtures:
    """Draw random polymer composites.

    Parameters
    ----------
    constituents : list of str
        Pool of shielding constituents (names in the constituent registry).
    n_constituents : int | list[int] | dict[int, float]
        Mixture size: ``2``; ``[1, 2, 3]`` (equally often); ``{2: 0.7, 3: 0.3}``.
    polymer_fractions : list of float
        Polymer weight fraction, drawn uniformly from this list.
    splitter : object with ``split(names, total, rng)``
        How the non-polymer mass is divided. Default: flat :class:`Dirichlet`.
    combinations : ``"balanced"`` | ``"uniform"``
        ``"balanced"`` spreads draws as evenly as possible over sizes and over
        the eligible combinations; ``"uniform"`` draws each mixture independently.
    require_any : list of str
        Only combinations containing at least one of these are eligible.
    decimals : int | None
        Round fractions to this many decimals, preserving the total and
        keeping every share non-zero. ``2`` matches the folder names (and the
        paper's data); ``None`` keeps full precision.
    """

    def __init__(self, constituents, *, n_constituents=2,
                 polymer_fractions=(0.20, 0.30, 0.40), polymer=POLYMER_FORMULA,
                 splitter=None, combinations="balanced", require_any=(), decimals=2):
        self.constituents = sorted(set(constituents))
        self.size_weights = _normalise_weights(n_constituents)
        self.polymer_fractions = [float(p) for p in polymer_fractions]
        self.polymer = polymer
        self.splitter = splitter or Dirichlet()
        if combinations not in ("balanced", "uniform"):
            raise ValueError("combinations must be 'balanced' or 'uniform'")
        self.combinations = combinations
        self.require_any = set(require_any)
        self.decimals = decimals

        if not self.polymer_fractions or any(not 0 <= p < 1 for p in self.polymer_fractions):
            raise ValueError("polymer_fractions must be a non-empty list of values in [0, 1)")
        unknown = self.require_any - set(self.constituents)
        if unknown:
            raise ValueError(f"require_any lists constituents not in the pool: {sorted(unknown)}")
        self._combos = {}
        for k in self.size_weights:
            if not 1 <= k <= len(self.constituents):
                raise ValueError(f"Cannot pick {k} constituents from a pool of "
                                 f"{len(self.constituents)}")
            combos = [c for c in itertools.combinations(self.constituents, k)
                      if not self.require_any or self.require_any.intersection(c)]
            if not combos:
                raise ValueError(f"No {k}-constituent combination contains any of "
                                 f"{sorted(self.require_any)}")
            self._combos[k] = combos

    def _schedule(self, count, rng):
        if self.combinations == "uniform":
            sizes = list(self.size_weights)
            ks = rng.choices(sizes, weights=[self.size_weights[k] for k in sizes], k=count)
            return [rng.choice(self._combos[k]) for k in ks]
        schedule = []
        for k, n_k in sorted(_allocate(count, self.size_weights).items()):
            combos = self._combos[k][:]
            rng.shuffle(combos)
            schedule += [combos[i % len(combos)] for i in range(n_k)]
        rng.shuffle(schedule)
        return schedule

    def mixture(self, names, rng):
        """One mixture of exactly these constituents."""
        poly = rng.choice(self.polymer_fractions)
        fracs = self.splitter.split(list(names), 1.0 - poly, rng)
        if self.decimals is not None:
            fracs = round_preserving_sum(fracs, self.decimals)
        return Mixture(dict(zip(names, fracs)), poly, self.polymer)

    def draw(self, count, rng):
        """``count`` mixtures using the random generator ``rng``."""
        return [self.mixture(names, rng) for names in self._schedule(count, rng)]

    def sample(self, count, seed=None):
        """``count`` mixtures from a fresh generator seeded with ``seed``."""
        return self.draw(count, random.Random(seed))


# ---------------------------------------------------------------------------
# Samplers
# ---------------------------------------------------------------------------
def thickness_range(start_mm, stop_mm, step_mm):
    """Inclusive, rounding-safe range: ``thickness_range(0.15, 1.5, 0.15)``."""
    n = int(round((stop_mm - start_mm) / step_mm))
    return [round(start_mm + i * step_mm, 6) for i in range(n + 1)]


class RandomSampler(Sampler):
    """``n`` random mixtures per beam, each at one randomly chosen thickness.

    ``n`` is a count per beam (with ``beams``) or ``{beam: count}``. Pass
    either the :class:`RandomMixtures` keyword arguments directly
    (``constituents=...``, ``n_constituents=...``, ...) or a ready-made
    ``mixtures=`` object with a ``draw(count, rng)`` method.

    ``unique`` avoids repeats by redrawing: ``"spec"`` (default) never emits
    the same simulation twice; ``"mixture"`` never repeats a composition at
    all, not even at another thickness or beam (the paper's rule); ``None``
    allows repeats. If the design space is too small, a warning says how
    many repeats remain.

        RandomSampler(n=300, beams=["RQR6_080", "RQR8_100"],
                      constituents=["Bi2O3", "BaSO4", "W", "Sn"],
                      n_constituents=[1, 2, 3],
                      polymer_fractions=[0.2, 0.3, 0.4],
                      thicknesses_mm=[0.15, 0.3, 0.45, 0.6, 0.75],
                      seed=0)
    """

    def __init__(self, n, beams=None, constituents=None, *, thicknesses_mm,
                 mixtures=None, seed=None, penelope_seed=(1, 1), unique="spec",
                 max_redraws=100, **mixture_kwargs):
        if unique not in ("spec", "mixture", None):
            raise ValueError("unique must be 'spec', 'mixture' or None")
        if hasattr(n, "items"):
            if beams:
                raise ValueError("Give either beams=[...] with n=<int>, or n={beam: count}")
            self.counts = dict(n)
        else:
            if not beams:
                raise ValueError("Give beams=[...] or n={beam: count}")
            self.counts = {b: int(n) for b in beams}
        if mixtures is None:
            if constituents is None:
                raise ValueError("Give constituents=[...] or mixtures=RandomMixtures(...)")
            mixtures = RandomMixtures(constituents, **mixture_kwargs)
        elif mixture_kwargs or constituents is not None:
            raise ValueError("Pass either mixtures= or RandomMixtures keyword arguments, not both")
        self.mixtures = mixtures
        self.thicknesses_mm = [float(t) for t in thicknesses_mm]
        if not self.thicknesses_mm:
            raise ValueError("thicknesses_mm is empty")
        self.seed = seed
        self.penelope_seed = penelope_seed
        self.unique = unique
        self.max_redraws = max_redraws

    def _key(self, beam, mix, t):
        return mix if self.unique == "mixture" else (beam, mix, t)

    def __iter__(self):
        rng = random.Random(self.seed)
        seen, repeats = set(), 0
        for beam, count in self.counts.items():
            for mix in self.mixtures.draw(count, rng):
                t = rng.choice(self.thicknesses_mm)
                if self.unique:
                    for _ in range(self.max_redraws):
                        if self._key(beam, mix, t) not in seen:
                            break
                        if hasattr(self.mixtures, "mixture"):  # same constituents, new split
                            mix = self.mixtures.mixture(mix.names, rng)
                        t = rng.choice(self.thicknesses_mm)
                    key = self._key(beam, mix, t)
                    repeats += key in seen
                    seen.add(key)
                yield SimulationSpec(beam, mix, t, seed=self.penelope_seed,
                                     tags={"design": "random"})
        if repeats:
            log.warning("RandomSampler: %d repeated %s(s) remain after %d redraws each; the "
                        "design space is too small for unique=%r", repeats, self.unique,
                        self.max_redraws, self.unique)

    def __len__(self):
        return sum(self.counts.values())


class GridSampler(Sampler):
    """Every combination of ``beams`` x ``mixtures`` x ``thicknesses_mm``."""

    def __init__(self, beams, mixtures, thicknesses_mm, *, penelope_seed=(1, 1),
                 design="grid"):
        self.beams = list(beams)
        self.mixtures = [mixtures] if isinstance(mixtures, Mixture) else list(mixtures)
        self.thicknesses_mm = [float(t) for t in thicknesses_mm]
        self.penelope_seed = penelope_seed
        self.design = design

    def __iter__(self):
        for beam, mix, t in itertools.product(self.beams, self.mixtures, self.thicknesses_mm):
            yield SimulationSpec(beam, mix, t, seed=self.penelope_seed,
                                 tags={"design": self.design})

    def __len__(self):
        return len(self.beams) * len(self.mixtures) * len(self.thicknesses_mm)


class ThicknessSweep(GridSampler):
    """Each mixture at every thickness in the list, for each beam.

        ThicknessSweep(Mixture({"Bi2O3": 0.4, "BaSO4": 0.2}, 0.4),
                       thickness_range(0.15, 1.5, 0.15), beams=["RQR6_080"])

    ``mixtures`` may be one :class:`Mixture`, a list, or random ones:
    ``RandomMixtures([...]).sample(6, seed=1)``.
    """

    def __init__(self, mixtures, thicknesses_mm, beams, **kwargs):
        kwargs.setdefault("design", "sweep")
        super().__init__(beams, mixtures, thicknesses_mm, **kwargs)


def wide_csv_row(row, max_molecules=10, polymer=POLYMER_FORMULA):
    """Parse one row of a wide-format CSV (the layout :class:`WideColumns` writes).

    Columns: ``molecule_i`` / ``ratio_i`` (percent, ``-1`` = unused),
    ``polymer_ratio``, ``geo_number_mm`` (hundredths of a mm), ``spc``
    (tube kV like ``80`` or a beam key).
    """
    fractions = {}
    for i in range(1, max_molecules + 1):
        name = (row.get(f"molecule_{i}") or "-1").strip()
        if name == "-1":
            continue
        ratio = float(row[f"ratio_{i}"])
        if ratio > 0:
            fractions[name] = ratio / 100.0
    if not fractions:
        return None
    spc = row["spc"].strip()
    beam = f"{float(spc):g}" if spc.replace(".", "", 1).isdigit() else spc
    return SimulationSpec(
        beam=beam,
        mixture=Mixture(fractions, float(row["polymer_ratio"]), polymer),
        thickness_mm=float(row["geo_number_mm"]) / 100.0,
        tags={"design": "csv"},
    )


class CSVSampler(Sampler):
    """Specs from a CSV file, one row at a time.

    ``parse_row(row_dict)`` returns a spec, an iterable of specs, or ``None``
    to skip the row. The default reads the wide-format CSV layout (see
    :func:`wide_csv_row`).
    """

    def __init__(self, path, parse_row=wide_csv_row, **reader_kwargs):
        self.path = path
        self.parse_row = parse_row
        self.reader_kwargs = reader_kwargs

    def __iter__(self):
        with open(self.path, newline="") as f:
            for row in csv.DictReader(f, **self.reader_kwargs):
                out = self.parse_row(row)
                if out is None:
                    continue
                if isinstance(out, SimulationSpec):
                    yield out
                else:
                    yield from out
