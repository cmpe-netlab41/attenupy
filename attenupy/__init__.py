"""attenupy -- design, build, run and harvest penEasy/PENELOPE shielding studies.

You supply two executables compiled for your machine (``penEasy.x`` and
PENELOPE's ``material.exe``); attenupy generates everything else: folders,
``penEasy.in``, the ``.geo`` geometry, the ``.spc`` spectrum, ``mat.mat`` and
``Aire.mat``.

    import attenupy as atp

    design = atp.RandomSampler(
        n=300, beams=["RQR6_080", "RQR8_100"],
        constituents=["Bi2O3", "BaSO4", "W", "Sn"], n_constituents=[1, 2, 3],
        polymer_fractions=[0.2, 0.3, 0.4],
        thicknesses_mm=atp.thickness_range(0.15, 1.5, 0.15), seed=0,
    )
    study = atp.Study("runs/demo", toolchain=atp.Toolchain.from_env())
    study.prepare(design)
    study.run(max_concurrent=10)
    study.collect(output="results.csv")

Extension points (duck-typed; subclass nothing unless you want to):

=====================  ===============================================  =======================
what                   interface                                        built in
=====================  ===============================================  =======================
design                 iterable of SimulationSpec (or ``Sampler``)      RandomSampler, GridSampler,
                                                                        ThicknessSweep, CSVSampler
composition            ``draw(count, rng) -> [Mixture]``                RandomMixtures
mass split             ``split(names, total, rng) -> [float]``          Dirichlet, EqualSplit
geometry               ``filename(spec)``, ``write(spec, path)``        SlabGeometry, FileGeometry
spectrum               ``write(path)``                                  FileSpectrum, TabulatedSpectrum,
                                                                        SpekpySpectrum
materials              ``build(mixture, path, reg)``, ``build_air(p)``  PenelopeMaterials
folder names           ``name(spec)``                                   CompactNamer, DescriptiveNamer
runner                 ``run(dirs)``                                    LocalRunner
result columns         ``extractor(sim, ctx) -> dict``                  WideColumns, SpecColumns,
                                                                        material_properties, EnergyRatio
=====================  ===============================================  =======================
"""

__version__ = "1.0.0"

from .beams import DEFAULT_BEAMS, Beam, BeamRegistry, register_beam
from .chemistry import (build_mixture, element_weight_fractions, mixture_density,
                        molar_mass, parse_formula)
from .constituents import (DEFAULT_CONSTITUENTS, Constituent, ConstituentRegistry,
                           register_constituent)
from .data import atomic_numbers, atomic_weights, densities
from .design import load_design, register_sampler
from .geometry import FileGeometry, SlabGeometry
from .materials import MaterialToolError, PenelopeMaterials, PenelopeMaterialTool
from .naming import (CompactNamer, DescriptiveNamer, mm_to_token,
                     parse_folder_name, token_to_mm)
from .results import (DEFAULT_EXTRACTORS, CollectContext, EnergyRatio, Simulation,
                      SpecColumns, WideColumns, collect, find_simulations, load_spec,
                      material_properties, mass_density, mean_excitation,
                      read_energy_deposition, write_csv)
from .runner import LocalRunner, RunSummary, sim_state
from .sampling import (Chain, CSVSampler, Dirichlet, EqualSplit, GridSampler,
                       RandomMixtures, RandomSampler, Sampler, ThicknessSweep,
                       round_preserving_sum, thickness_range)
from .spec import Mixture, SimulationSpec
from .spectra import (FileSpectrum, SpekpySpectrum, TabulatedSpectrum, bundled_spectra,
                      bundled_spectrum, check_bigs_spectrum, read_spectrum,
                      write_bigs_spectrum)
from .study import Study
from .template import PeneasyDeck
from .toolchain import Toolchain, ToolchainError


def log_to_console(level="INFO"):
    """Show attenupy's progress messages (it logs to the ``attenupy`` logger)."""
    import logging
    logger = logging.getLogger("attenupy")
    if not any(getattr(h, "_attenupy", False) for h in logger.handlers):
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(message)s"))
        handler._attenupy = True
        logger.addHandler(handler)
    logger.setLevel(level)
    logger.propagate = False


import types as _types

__all__ = sorted(n for n, v in globals().items()
                 if not n.startswith("_") and not isinstance(v, _types.ModuleType))
