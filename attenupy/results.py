"""Harvest finished simulations into rows (and CSV).

Collection is a list of *extractors*. An extractor is any callable::

    extractor(sim, ctx) -> dict        # column name -> value

where ``sim`` is a :class:`Simulation` (``.path``, ``.spec``, ``.beam``) and
``ctx`` a :class:`CollectContext` (reference energies, beam registry). Rows
are the merged dicts, in extractor order. For example, a custom column::

    def detector_uncertainty(sim, ctx):
        tally = read_energy_deposition(sim.file("tallyEnergyDeposition.dat"))
        return {"det_unc_2sigma": tally[2][1] if tally else None}

    study.collect(extractors=DEFAULT_EXTRACTORS + [detector_uncertainty])

:data:`DEFAULT_EXTRACTORS` produce a wide CSV row per simulation
(``molecule_*``, ``ratio_*``, ``polymer``, ``polymer_ratio``,
``geo_number_mm``, ``spc``, ``mass_den``, ``mean_excit``, ``energy_ratio``).
"""

import csv
import logging
import os
import re
from dataclasses import dataclass

from .beams import DEFAULT_BEAMS
from .constituents import DEFAULT_CONSTITUENTS
from .naming import CompactNamer
from .spec import SPEC_FILENAME, SimulationSpec

log = logging.getLogger(__name__)

DONE_MARKER = "Have a nice day."
ENERGY_TALLY = "tallyEnergyDeposition.dat"
DETECTOR_MATERIAL = 2

_MASS_DENSITY_RE = re.compile(r"Mass density = ([\d.Ee+-]+) g/cm\*\*3")
_MEAN_EXC_RE = re.compile(r"Mean excitation energy = ([\d.Ee+-]+) eV")


# ---------------------------------------------------------------------------
# File readers
# ---------------------------------------------------------------------------
def _scan(path, regex):
    try:
        with open(path, errors="replace") as f:
            for line in f:
                m = regex.search(line)
                if m:
                    return float(m.group(1))
    except OSError:
        pass
    return None


def mass_density(mat_file):
    """Mass density (g/cm**3) reported inside a PENELOPE ``.mat`` file."""
    return _scan(mat_file, _MASS_DENSITY_RE)


def mean_excitation(mat_file):
    """Mean excitation energy (eV) reported inside a PENELOPE ``.mat`` file."""
    return _scan(mat_file, _MEAN_EXC_RE)


def read_energy_deposition(tally_file):
    """``{material: (eV_per_history, two_sigma)}`` from penEasy's EDP tally."""
    out = {}
    try:
        with open(tally_file) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 3 and not line.lstrip().startswith("#"):
                    try:
                        out[int(parts[0])] = (float(parts[1]), float(parts[2]))
                    except ValueError:
                        continue
    except OSError:
        return None
    return out or None


def is_finished(sim_dir):
    """True when penEasy ran to completion (its output ends with the sign-off)."""
    path = os.path.join(sim_dir, "penEasy.out")
    try:
        with open(path, "rb") as f:
            f.seek(max(0, os.path.getsize(path) - 2048))
            return DONE_MARKER.encode() in f.read()
    except OSError:
        return False


# ---------------------------------------------------------------------------
# Simulations on disk
# ---------------------------------------------------------------------------
def load_spec(sim_dir):
    """The spec of a simulation directory: ``spec.json``, or parsed from its
    folder name. Raises ``ValueError`` if neither works."""
    if os.path.isfile(os.path.join(sim_dir, SPEC_FILENAME)):
        return SimulationSpec.load(sim_dir)
    return CompactNamer.parse(os.path.basename(os.path.normpath(sim_dir)))


@dataclass
class Simulation:
    path: str
    spec: SimulationSpec
    beam: object = None      # resolved Beam, if the registry knows it

    def file(self, name):
        return os.path.join(self.path, name)

    @property
    def finished(self):
        return is_finished(self.path)


def find_simulations(root, beams=None, include_references=False):
    """All simulation directories under ``root`` (any depth, sorted).

    A directory counts when it holds a ``penEasy.in``. Directories whose spec
    cannot be determined are skipped with a message.
    """
    beams = beams or DEFAULT_BEAMS
    sims = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(d for d in dirnames if d != ".attenupy")
        if "penEasy.in" not in filenames:
            continue
        dirnames[:] = []  # simulations don't nest
        try:
            spec = load_spec(dirpath)
        except (ValueError, KeyError) as e:
            log.warning("skipping %s: %s", dirpath, e)
            continue
        if spec.is_reference and not include_references:
            continue
        beam = beams.resolve(spec.beam) if spec.beam in beams else None
        sims.append(Simulation(dirpath, spec, beam))
    return sims


# ---------------------------------------------------------------------------
# Context + extractors
# ---------------------------------------------------------------------------
class CollectContext:
    """Shared state for extractors: where I0 comes from, warnings shown once."""

    def __init__(self, beams=None, measured_references=None, constituents=None,
                 density_tolerance=0.01):
        self.beams = beams or DEFAULT_BEAMS
        self.constituents = constituents or DEFAULT_CONSTITUENTS
        self.measured = dict(measured_references or {})
        self.density_tolerance = density_tolerance
        self.density_mismatches = []
        self._warned = set()

    def reference_energy(self, sim):
        """I0 for ``sim``'s beam: a measured air run if the study has one,
        otherwise the beam's configured ``reference_energy``."""
        key = sim.beam.key if sim.beam else sim.spec.beam
        if key in self.measured:
            return self.measured[key]
        if sim.beam is not None and sim.beam.reference_energy is not None:
            return sim.beam.reference_energy
        self.warn_once(key, f"No reference energy (I0) for beam {key!r}: energy_ratio "
                            "left empty. Run study.prepare_references() and simulate it.")
        return None

    def warn_once(self, key, message):
        if key not in self._warned:
            self._warned.add(key)
            log.warning("%s", message)


class WideColumns:
    """``molecule_1..N``, ``ratio_1..N`` (percent), ``polymer``,
    ``polymer_ratio``, ``geo_number_mm`` (hundredths of a mm), ``spc`` (kV)."""

    def __init__(self, max_molecules=10):
        self.max_molecules = max_molecules

    def __call__(self, sim, ctx):
        mix = sim.spec.mixture
        if mix.n_constituents > self.max_molecules:
            raise ValueError(f"{sim.path}: {mix.n_constituents} constituents > "
                             f"max_molecules={self.max_molecules}")
        pad = self.max_molecules - mix.n_constituents
        names = mix.names + ["-1"] * pad
        ratios = [f * 100.0 for _, f in mix.fractions] + [-1.0] * pad
        row = {f"molecule_{i + 1}": n for i, n in enumerate(names)}
        row.update({f"ratio_{i + 1}": r for i, r in enumerate(ratios)})
        row.update({
            "polymer": mix.polymer,
            "polymer_ratio": mix.polymer_fraction,
            "geo_number_mm": int(round(sim.spec.thickness_mm * 100)),
            "spc": int(sim.beam.kvp) if sim.beam else None,
        })
        return row


class SpecColumns:
    """Tidy description: beam, kvp, thickness_mm, polymer, polymer_fraction,
    n_constituents and one ``w_<name>`` column per constituent."""

    def __call__(self, sim, ctx):
        mix = sim.spec.mixture
        row = {
            "beam": sim.spec.beam,
            "kvp": sim.beam.kvp if sim.beam else None,
            "thickness_mm": sim.spec.thickness_mm,
            "polymer": mix.polymer,
            "polymer_fraction": mix.polymer_fraction,
            "n_constituents": mix.n_constituents,
        }
        row.update({f"w_{n}": f for n, f in mix.fractions})
        return row


def material_properties(sim, ctx):
    """``mass_den`` and ``mean_excit`` from the slab's ``mat.mat``.

    Also cross-checks the density against the spec: a mismatch means the
    ``mat.mat`` in the folder is not the material the spec describes (e.g. a
    stale file), and is recorded in ``ctx.density_mismatches``.
    """
    mat = sim.file("mat.mat")
    rho = mass_density(mat)
    if rho and sim.spec.mixture is not None:
        try:
            expected = sim.spec.mixture.density(ctx.constituents)
        except KeyError:
            expected = None
        if expected and abs(rho - expected) / expected > ctx.density_tolerance:
            ctx.density_mismatches.append((sim.path, rho, expected))
    return {"mass_den": rho, "mean_excit": mean_excitation(mat)}


class EnergyRatio:
    """``energy_ratio`` = detector energy / I0 (the transmission y = I/I0)."""

    def __init__(self, detector_material=DETECTOR_MATERIAL, column="energy_ratio"):
        self.detector_material = detector_material
        self.column = column

    def __call__(self, sim, ctx):
        tally = read_energy_deposition(sim.file(ENERGY_TALLY))
        i0 = ctx.reference_energy(sim)
        value = None
        if tally and self.detector_material in tally and i0:
            value = tally[self.detector_material][0] / i0
        return {self.column: value}


def path_column(sim, ctx):
    return {"path": sim.path}


DEFAULT_EXTRACTORS = [WideColumns(), material_properties, EnergyRatio()]


# ---------------------------------------------------------------------------
# Collection
# ---------------------------------------------------------------------------
def collect(simulations, extractors=None, ctx=None, *, only_finished=True):
    """Apply ``extractors`` to each :class:`Simulation`; return a list of rows."""
    extractors = DEFAULT_EXTRACTORS if extractors is None else extractors
    ctx = ctx or CollectContext()
    rows, skipped = [], 0
    for sim in simulations:
        if only_finished and not sim.finished:
            skipped += 1
            continue
        row = {}
        for extract in extractors:
            row.update(extract(sim, ctx))
        rows.append(row)
    log.info("Collected %d simulations%s.", len(rows),
             f" ({skipped} unfinished skipped)" if skipped else "")
    if ctx.density_mismatches:
        examples = "".join(f"\n  {p}: mat.mat {rho:.4f} g/cm3, spec {exp:.4f}"
                           for p, rho, exp in ctx.density_mismatches[:5])
        log.warning("%d mat.mat file(s) disagree with their spec by more than %.0f%% in "
                    "density, e.g.:%s", len(ctx.density_mismatches),
                    ctx.density_tolerance * 100, examples)
    return rows


def write_csv(rows, path, columns=None):
    """Write rows; columns default to the union of keys in first-seen order."""
    if columns is None:
        columns = list(dict.fromkeys(k for row in rows for k in row))
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in columns})
    return path
