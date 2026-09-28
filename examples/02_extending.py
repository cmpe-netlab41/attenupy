"""How to plug your own pieces into attenupy.

Every stage is duck-typed: pass any object with the right method(s). This
script only writes input files (no executables needed):

    pip install spekpy   # optional, used by section 1 below
    python examples/02_extending.py
"""

import importlib.util
import math
import random
import tempfile

import attenupy as atp

HAVE_SPEKPY = importlib.util.find_spec("spekpy") is not None


# --- 1. A new constituent and a new beam at a custom energy -----------------
atp.register_constituent("Gd2O3", density=7.41)
atp.register_constituent("PU", density=1.20, formula="C3H8N2O")  # label != formula


class ExponentialToySpectrum:
    """A spectrum source: any object with ``write(path)``.

    Fallback for when SpekPy isn't installed; prefer ``atp.SpekpySpectrum``
    for anything real.
    """

    def __init__(self, kvp):
        self.kvp = kvp

    def write(self, path):
        energies = [e * 1e3 for e in range(10, int(self.kvp) + 1)]
        weights = [(self.kvp - e / 1e3) * math.exp(-e / 25e3) for e in energies]
        atp.write_bigs_spectrum(path, energies, weights, header=f"toy {self.kvp} kVp")


if HAVE_SPEKPY:
    spectrum_90kv = atp.SpekpySpectrum(kvp=90, th=12, filters=[("Al", 3.0), ("Be", 3.0)])
else:
    spectrum_90kv = ExponentialToySpectrum(90)

atp.register_beam("W090", kvp=90, spectrum=spectrum_90kv)
demo_beams = ["W090", "RQR6_080"]


# --- 2. A custom mass splitter: split(names, total, rng) -> list ------------
class HeavyFirst:
    """Give the heaviest-density constituent 60-80 % of the non-polymer mass."""

    def split(self, names, total, rng):
        if len(names) == 1:
            return [total]
        dens = [atp.DEFAULT_CONSTITUENTS[n].density for n in names]
        heavy = dens.index(max(dens))
        share = rng.uniform(0.6, 0.8)
        rest = (1 - share) / (len(names) - 1)
        return [total * (share if i == heavy else rest) for i in range(len(names))]


# --- 3. A custom sampler: subclass Sampler and implement __iter__ -----------
class LatinThickness(atp.Sampler):
    """Stratified thicknesses: one draw from each of ``n`` equal-width bins."""

    def __init__(self, mixture, beam, lo, hi, n, seed=0):
        self.mixture, self.beam, self.lo, self.hi, self.n, self.seed = mixture, beam, lo, hi, n, seed

    def __iter__(self):
        rng = random.Random(self.seed)
        width = (self.hi - self.lo) / self.n
        for i in range(self.n):
            t = round(self.lo + (i + rng.random()) * width, 2)
            yield atp.SimulationSpec(self.beam, self.mixture, t, tags={"design": "latin"})


# --- 4. A custom result column: extractor(sim, ctx) -> dict -----------------
def detector_uncertainty(sim, ctx):
    tally = atp.read_energy_deposition(sim.file("tallyEnergyDeposition.dat"))
    return {"det_2sigma": tally[2][1] if tally else None}


# --- Put it together ---------------------------------------------------------
design = (
    atp.RandomSampler(n=5, beams=demo_beams, constituents=["Gd2O3", "W", "Bi2O3"],
                     n_constituents=2, splitter=HeavyFirst(),
                     polymer_fractions=[0.25, 0.35],        # not multiples of 0.1 ...
                     thicknesses_mm=[0.5, 1.0], seed=1)
    + LatinThickness(atp.Mixture({"Gd2O3": 0.5}, 0.5), "RQR6_080", 0.1, 2.0, n=4)
)

root = tempfile.mkdtemp(prefix="attenupy-demo-")
study = atp.Study(
    root,
    namer=atp.DescriptiveNamer(),                      # ... so use unambiguous names
    geometry=atp.SlabGeometry(slab_half_width_cm=5.0),  # 10x10 cm slab
    deck=atp.PeneasyDeck(num_histories="1.0e8"),
)
paths = study.prepare(design, build_materials=False)  # inputs only; no material.exe here
for p in paths:
    print(p[len(root) + 1:])
print("status:", dict(study.status()))

# After running, collect with the extra column:
#   study.collect(atp.DEFAULT_EXTRACTORS + [detector_uncertainty], output="results.csv")
