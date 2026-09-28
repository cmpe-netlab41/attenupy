"""What to simulate: :class:`Mixture` (a composition) and :class:`SimulationSpec`.

A spec is a plain, immutable value. Samplers produce specs, a
:class:`~attenupy.study.Study` turns them into directories, and every
directory stores its spec as ``spec.json`` so later stages (materials, run,
collect) never have to reverse-engineer parameters from the folder name.
"""

import json
import math
import os
from dataclasses import dataclass, field, replace

from .chemistry import element_weight_fractions
from .constituents import DEFAULT_CONSTITUENTS
from .data import POLYMER_FORMULA

SPEC_FILENAME = "spec.json"

# Tolerance on |sum(fractions) - 1|. Loose on purpose: folder names store
# fractions rounded to 2 decimals, so their sums can be off by ~0.01.
FRACTION_SUM_TOL = 0.02


# RANECU (L'Ecuyer 1988) seed ranges. penEasy gives 0,0 (read a seeds file) and
# -1,-1 (restart from a dump) special meanings, so those are rejected too.
_SEED_MAX = (2147483562, 2147483398)


def _as_seed(seed):
    s1, s2 = (seed, seed) if isinstance(seed, int) else seed
    s1, s2 = int(s1), int(s2)
    if not (1 <= s1 <= _SEED_MAX[0] and 1 <= s2 <= _SEED_MAX[1]):
        raise ValueError(
            f"penEasy seeds must be in [1, {_SEED_MAX[0]}] and [1, {_SEED_MAX[1]}], "
            f"got ({s1}, {s2}); 0,0 and -1,-1 mean 'seeds file' and 'restart' to penEasy"
        )
    return (s1, s2)


@dataclass(frozen=True)
class Mixture:
    """A polymer-matrix composite: shielding constituents + polymer.

    ``fractions`` maps constituent name -> weight fraction *of the whole
    mixture*; together with ``polymer_fraction`` they sum to 1::

        Mixture({"Bi2O3": 0.40, "BaSO4": 0.20}, polymer_fraction=0.40)

    Stored sorted by name, so two mixtures with the same content compare equal
    (and hash equal, which is what the material cache relies on).
    """

    fractions: tuple
    polymer_fraction: float = 0.0
    polymer: str = POLYMER_FORMULA

    def __post_init__(self):
        items = self.fractions.items() if hasattr(self.fractions, "items") else self.fractions
        items = tuple(sorted((str(n), float(f)) for n, f in items))
        object.__setattr__(self, "fractions", items)
        object.__setattr__(self, "polymer_fraction", float(self.polymer_fraction))

        names = [n for n, _ in items]
        if not items:
            raise ValueError("A mixture needs at least one non-polymer constituent")
        if len(set(names)) != len(names):
            raise ValueError(f"Duplicate constituents in {names}")
        if self.polymer in names:
            raise ValueError(f"The polymer {self.polymer!r} is also listed as a constituent")
        if any(f <= 0 for _, f in items):
            raise ValueError(f"All fractions must be > 0, got {dict(items)}")
        if not 0.0 <= self.polymer_fraction < 1.0:
            raise ValueError(f"polymer_fraction must be in [0, 1), got {self.polymer_fraction}")
        total = self.polymer_fraction + sum(f for _, f in items)
        if abs(total - 1.0) > FRACTION_SUM_TOL:
            raise ValueError(
                f"Fractions must sum to 1 (polymer included), got {total:.4f} for "
                f"{dict(items)} + {self.polymer}={self.polymer_fraction}"
            )

    # -- views -------------------------------------------------------------
    def as_dict(self):
        """Non-polymer constituents as ``{name: fraction}``."""
        return dict(self.fractions)

    @property
    def names(self):
        return [n for n, _ in self.fractions]

    @property
    def n_constituents(self):
        return len(self.fractions)

    def components(self):
        """``[(name, weight), ...]`` with the polymer first (when present)."""
        comps = [(self.polymer, self.polymer_fraction)] if self.polymer_fraction > 0 else []
        return comps + list(self.fractions)

    @property
    def label(self):
        """``"BaSO4_Bi2O3_0.20_0.40"``: names then 2-decimal fractions."""
        return "_".join(self.names + [f"{f:.2f}" for _, f in self.fractions])

    # -- physics -----------------------------------------------------------
    def element_fractions(self, registry=None):
        """Element weight fractions ``{symbol: w}`` of the whole mixture."""
        registry = registry or DEFAULT_CONSTITUENTS
        return element_weight_fractions(
            (registry[name].formula, w) for name, w in self.components()
        )

    def density(self, registry=None):
        """Bulk density (g/cm**3) by the additivity-of-volumes rule."""
        registry = registry or DEFAULT_CONSTITUENTS
        return 1.0 / sum(w / registry[name].density for name, w in self.components())

    # -- (de)serialisation ------------------------------------------------
    def to_dict(self):
        return {
            "fractions": self.as_dict(),
            "polymer": self.polymer,
            "polymer_fraction": self.polymer_fraction,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(d["fractions"], d.get("polymer_fraction", 0.0),
                   d.get("polymer", POLYMER_FORMULA))


@dataclass(frozen=True)
class SimulationSpec:
    """One penEasy run: a beam, a slab of ``mixture`` ``thickness_mm`` thick.

    ``beam`` is a beam key or alias (``"RQR6_080"``, ``"80"``, ``"80kV"``); the
    study resolves it against its beam registry. ``kind="reference"`` marks an
    air run (the slab replaced by air) used to measure the unattenuated
    detector signal I0; it has no mixture.
    """

    beam: str
    mixture: Mixture = None
    thickness_mm: float = 0.0
    seed: tuple = (1, 1)
    kind: str = "sample"
    tags: dict = field(default_factory=dict, compare=False, hash=False)

    def __post_init__(self):
        object.__setattr__(self, "seed", _as_seed(self.seed))
        object.__setattr__(self, "thickness_mm", float(self.thickness_mm))
        if self.kind not in ("sample", "reference"):
            raise ValueError(f"kind must be 'sample' or 'reference', got {self.kind!r}")
        if self.kind == "sample" and self.mixture is None:
            raise ValueError("A sample spec needs a mixture")
        if not math.isfinite(self.thickness_mm) or self.thickness_mm <= 0:
            raise ValueError(f"thickness_mm must be > 0, got {self.thickness_mm}")

    @classmethod
    def reference(cls, beam, thickness_mm=0.30, seed=(1, 1)):
        """The air run for ``beam`` (the slab filled with air)."""
        return cls(beam=beam, mixture=None, thickness_mm=thickness_mm, seed=seed,
                   kind="reference")

    @property
    def is_reference(self):
        return self.kind == "reference"

    def replace(self, **changes):
        """A copy with some fields changed, e.g. ``spec.replace(thickness_mm=1.5)``."""
        return replace(self, **changes)

    # -- (de)serialisation ------------------------------------------------
    def to_dict(self):
        d = {
            "kind": self.kind,
            "beam": self.beam,
            "thickness_mm": self.thickness_mm,
            "seed": list(self.seed),
        }
        if self.mixture is not None:
            d.update(self.mixture.to_dict())
        if self.tags:
            d["tags"] = dict(self.tags)
        return d

    @classmethod
    def from_dict(cls, d):
        mixture = Mixture.from_dict(d) if "fractions" in d else None
        return cls(
            beam=d["beam"],
            mixture=mixture,
            thickness_mm=d["thickness_mm"],
            seed=tuple(d.get("seed", (1, 1))),
            kind=d.get("kind", "sample"),
            tags=dict(d.get("tags", {})),
        )

    def save(self, directory):
        path = os.path.join(directory, SPEC_FILENAME)
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2, sort_keys=True)
            f.write("\n")
        return path

    @classmethod
    def load(cls, directory):
        with open(os.path.join(directory, SPEC_FILENAME)) as f:
            return cls.from_dict(json.load(f))
