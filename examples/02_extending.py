"""How to plug your own pieces into attenupy.

attenupy is built around a few small, duck-typed interfaces: for each stage
of a study (spectrum, mass split, design, result column, ...) you can pass
*any* object that has the right method, instead of one of attenupy's
built-in ones. This script shows one custom implementation for each of those
extension points, with comments explaining the interface it satisfies. It is
meant to be read as much as run.

This script only writes input files (no penEasy.x or material.exe needed):

    pip install spekpy   # optional, used by section 1 below
    python examples/02_extending.py

If SpekPy is not installed, section 1 falls back to a hand-written toy
spectrum instead, so the script still runs end to end either way.
"""

import importlib.util
import math
import random
import tempfile

import attenupy as atp

# Check once, at import time, whether the optional SpekPy package is
# available, so section 1 below can pick a spectrum source without crashing.
HAVE_SPEKPY = importlib.util.find_spec("spekpy") is not None


# --- 1. A new constituent and a new beam at a custom energy -----------------
# register_constituent teaches attenupy about a shielding compound it
# doesn't already know: give it a name, a density (g/cm^3), and, if the name
# isn't itself a valid chemical formula, an explicit formula= for computing
# element weight fractions.
atp.register_constituent("Gd2O3", density=7.41)
atp.register_constituent("PU", density=1.20, formula="C3H8N2O")  # label != formula


class ExponentialToySpectrum:
    """A spectrum source: any object with ``write(path)``.

    This is the *interface* every spectrum source must implement: a
    ``write(path)`` method that writes a penEasy-format ``.spc`` file at
    ``path``. attenupy ships two ready-made spectrum sources
    (``FileSpectrum``, which copies an existing ``.spc`` file, and
    ``SpekpySpectrum``, which computes a real tube spectrum with SpekPy);
    this class is a third, hand-rolled one, used here only as a fallback
    for when SpekPy isn't installed. Prefer ``atp.SpekpySpectrum`` for
    anything real -- this one does not model an actual X-ray tube, it just
    produces *some* plausible-looking spectrum shape for the demo.
    """

    def __init__(self, kvp):
        self.kvp = kvp  # tube voltage in kV; also the maximum photon energy

    def write(self, path):
        # A made-up exponential falloff from 10 keV up to kvp keV -- good
        # enough to produce *a* spectrum for the demo, not a physically
        # accurate one.
        energies = [e * 1e3 for e in range(10, int(self.kvp) + 1)]  # eV
        weights = [(self.kvp - e / 1e3) * math.exp(-e / 25e3) for e in energies]
        atp.write_bigs_spectrum(path, energies, weights, header=f"toy {self.kvp} kVp")


# Pick a real SpekPy-computed spectrum when available, otherwise fall back
# to the toy one above -- either way we end up with one 90 kV beam, an
# energy none of attenupy's five built-in beams provide.
if HAVE_SPEKPY:
    spectrum_90kv = atp.SpekpySpectrum(kvp=90, th=12, filters=[("Al", 3.0), ("Be", 3.0)])
else:
    spectrum_90kv = ExponentialToySpectrum(90)

atp.register_beam("W090", kvp=90, spectrum=spectrum_90kv)
demo_beams = ["W090", "RQR6_080"]  # our new beam, plus one of the built-in ones


# --- 2. A custom mass splitter: split(names, total, rng) -> list ------------
class HeavyFirst:
    """Give the heaviest-density constituent 60-80 % of the non-polymer mass.

    A splitter decides how the non-polymer mass of a mixture is divided
    between its constituents. attenupy's default is ``Dirichlet()`` (a
    random, roughly even split); this one is a deliberately biased
    alternative that always favours whichever constituent is densest.

    ``names`` is the list of constituent names in this mixture, ``total`` is
    the non-polymer weight fraction to distribute among them (e.g. 0.75 if
    the polymer fraction is 0.25), and ``rng`` is a ``random.Random``
    instance to draw from. The returned list of fractions must sum to
    ``total`` and have one entry per name, in the same order as ``names``.
    """

    def split(self, names, total, rng):
        if len(names) == 1:
            return [total]
        dens = [atp.DEFAULT_CONSTITUENTS[n].density for n in names]
        heavy = dens.index(max(dens))           # index of the densest constituent
        share = rng.uniform(0.6, 0.8)           # its share of the total mass
        rest = (1 - share) / (len(names) - 1)   # the rest split evenly among the others
        return [total * (share if i == heavy else rest) for i in range(len(names))]


# --- 3. A custom sampler: subclass Sampler and implement __iter__ -----------
class LatinThickness(atp.Sampler):
    """Stratified thicknesses: one draw from each of ``n`` equal-width bins.

    A sampler is any iterable of ``SimulationSpec`` objects; subclassing
    ``atp.Sampler`` just adds the ``+`` operator for combining it with other
    samplers (see "Put it together" below). Instead of drawing thicknesses
    uniformly at random over the whole [lo, hi] range -- which can leave
    gaps or clusters purely by chance -- this splits the range into ``n``
    equal-width bins and draws one random thickness from each bin, a simple
    form of stratified sampling that spreads the points out more evenly.
    """

    def __init__(self, mixture, beam, lo, hi, n, seed=0):
        self.mixture, self.beam, self.lo, self.hi, self.n, self.seed = mixture, beam, lo, hi, n, seed

    def __iter__(self):
        rng = random.Random(self.seed)
        width = (self.hi - self.lo) / self.n
        for i in range(self.n):
            # One random point inside bin i: [lo + i*width, lo + (i+1)*width).
            t = round(self.lo + (i + rng.random()) * width, 2)
            yield atp.SimulationSpec(self.beam, self.mixture, t, tags={"design": "latin"})


# --- 4. A custom result column: extractor(sim, ctx) -> dict -----------------
def detector_uncertainty(sim, ctx):
    """Add a ``det_2sigma`` column with the detector tally's own uncertainty.

    An extractor is any callable ``extractor(sim, ctx) -> dict``; ``collect``
    runs a list of them over every finished simulation and merges their
    dicts into one output row per simulation. attenupy's built-in extractors
    (``WideColumns``, ``material_properties``, ``EnergyRatio``, ...) cover
    the usual columns; this shows how to add one of your own by reading a
    raw penEasy tally file directly.
    """
    tally = atp.read_energy_deposition(sim.file("tallyEnergyDeposition.dat"))
    return {"det_2sigma": tally[2][1] if tally else None}


# --- Put it together ---------------------------------------------------------
# Combine two different kinds of sampler with `+`: 5 random 2-constituent
# mixtures per beam (using our custom HeavyFirst splitter), plus 4
# stratified thicknesses of one fixed Gd2O3 mixture at 80 kV.
design = (
    atp.RandomSampler(n=5, beams=demo_beams, constituents=["Gd2O3", "W", "Bi2O3"],
                     n_constituents=2, splitter=HeavyFirst(),
                     polymer_fractions=[0.25, 0.35],        # not multiples of 0.1 ...
                     thicknesses_mm=[0.5, 1.0], seed=1)
    + LatinThickness(atp.Mixture({"Gd2O3": 0.5}, 0.5), "RQR6_080", 0.1, 2.0, n=4)
)

# A throwaway directory for this demo; use a real path in your own scripts.
root = tempfile.mkdtemp(prefix="attenupy-demo-")
study = atp.Study(
    root,
    # CompactNamer (the default) would refuse the 0.25/0.35 polymer
    # fractions above, since they aren't multiples of 0.1 and it would
    # mislabel them. DescriptiveNamer spells out every fraction explicitly
    # instead, so it accepts any fraction -- use it for designs like this one.
    namer=atp.DescriptiveNamer(),
    geometry=atp.SlabGeometry(slab_half_width_cm=5.0),  # a wider, 10x10 cm slab
    deck=atp.PeneasyDeck(num_histories="1.0e8"),         # fewer histories than the default
)

# build_materials=False writes only the text inputs (penEasy.in, .geo, .spc)
# and skips calling material.exe -- handy here since this demo has none.
paths = study.prepare(design, build_materials=False)
for p in paths:
    print(p[len(root) + 1:])   # each folder's path, relative to root

# --- 5. I0 for a brand-new beam -----------------------------------------
# W090 is a beam we just registered, so it has no baked-in reference energy
# (I0) -- only the five beams attenupy ships with do. Without one, collect()
# would warn and leave energy_ratio empty for every W090 simulation above.
# Add one air-reference simulation for it: same recipe as any other folder
# (prepare(), then run()), but with the slab filled with air instead of a
# real mixture, so its detector tally *is* I0 (see the README's "I0 and the
# transmission ratio"). It costs nothing to add for RQR6_080 too, but that
# built-in beam already has a reference energy, so it's left out here.
study.prepare_references(["W090"], build_materials=False)

print("status:", dict(study.status()))   # every folder, including the new air run, is "unprepared"

# After running penEasy.x for real (e.g. with build_materials=True and a
# configured Toolchain) -- including the air-reference folder just added --
# collect the usual columns plus our custom one. W090's energy_ratio will
# come from that reference run's measured I0 automatically:
#   study.collect(atp.DEFAULT_EXTRACTORS + [detector_uncertainty], output="results.csv")
