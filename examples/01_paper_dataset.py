"""Build the paper's dataset design end to end.

This script is a complete, runnable example: it defines a study design
(which simulations to run), creates the folders for every simulation, runs
penEasy, and collects the results into a CSV -- all in one pass. Read it top
to bottom as a tour of the four stages every attenupy study goes through:

    1. design   -- decide *what* to simulate (RandomSampler, ThicknessSweep, ...)
    2. prepare  -- turn that design into folders with all the input files
    3. run      -- execute penEasy.x for every folder
    4. collect  -- harvest the results into one CSV

You need your own penEasy.x and PENELOPE material.exe (attenupy does not
compile anything); point the environment variables at them before running:

    export ATTENUPY_PENEASY=/path/to/penEasy.x
    export ATTENUPY_MATERIAL=/path/to/pendbase/material.exe
    python examples/01_paper_dataset.py

This creates a `runs/paper/` directory with one subfolder per simulation and,
since the executables are configured, actually runs all of them -- for the
full 1200+ simulations below that can take a long time, so lower `n` if you
just want to see the mechanics work.
"""

import attenupy as atp

# Print progress messages (folders created, materials built, simulations
# launched, ...) to the terminal as the script runs. Without this the
# library stays silent unless something goes wrong.
atp.log_to_console()

# The twelve shielding compounds/elements used in the paper, mixed into a
# polyethylene (C2H4) matrix. These are just names attenupy already knows
# the density and chemical formula for (see atp.DEFAULT_CONSTITUENTS); if
# you wanted a constituent that isn't built in, you'd register it first with
# atp.register_constituent(name, density=...).
SHIELDERS = ["Bi2O3", "BaSO4", "BN", "B4C", "WO3", "Sb", "Sn", "Bi", "W", "Cu", "B", "Ba"]

# The candidate slab thicknesses, in millimetres: 0.15, 0.30, 0.45, ..., 1.50.
# thickness_range(start, stop, step) is just a shorthand for writing that
# list out by hand.
THICKNESSES = atp.thickness_range(0.15, 1.5, 0.15)  # mm

# --- Training set -----------------------------------------------------------
# For each of the 4 tube voltages ("beams" below), draw 300 random
# composites. For each one:
#   - pick 1, 2 or 3 shielding constituents out of SHIELDERS (n_constituents)
#   - pick a polymer weight fraction from {0.20, 0.30, 0.40}
#   - split the remaining (non-polymer) mass between the chosen constituents
#     with a flat Dirichlet distribution (the default mass splitter)
#   - pick one random thickness from THICKNESSES
# unique="mixture" makes sure the same composition is never generated twice,
# even at a different thickness or beam (the rule the paper's dataset uses).
# seed=2026 makes the random draw reproducible: running this script again
# with the same seed produces the exact same 1200 (4 beams x 300) composites.
training = atp.RandomSampler(
    n=300,
    beams=["N30_030", "RQR6_080", "RQR8_100", "x120kV"],  # 30, 80, 100, 120 kV
    constituents=SHIELDERS,
    n_constituents=[1, 2, 3],
    polymer_fractions=[0.20, 0.30, 0.40],
    thicknesses_mm=THICKNESSES,
    unique="mixture",
    seed=2026,
)

# --- Held-out curves ---------------------------------------------------------
# Two specific, hand-picked compositions, each simulated at *every*
# thickness in THICKNESSES (not just one random thickness, like above). This
# traces out a full transmission-vs-thickness curve for a fixed composition,
# useful for checking how well a model trained on `training` interpolates
# and extrapolates, on data it never saw during training.
held_out = atp.ThicknessSweep(
    [atp.Mixture({"Bi2O3": 0.40, "BaSO4": 0.20}, polymer_fraction=0.40),
     atp.Mixture({"Bi": 0.35, "W": 0.35}, polymer_fraction=0.30)],
    THICKNESSES,
    beams=["RQR6_080"],  # 80 kV only, for both compositions
)

# A Study ties a design to a directory on disk. toolchain=... tells it where
# penEasy.x and material.exe live; Toolchain.from_env() reads that from the
# ATTENUPY_PENEASY / ATTENUPY_MATERIAL environment variables set above.
study = atp.Study("runs/paper", toolchain=atp.Toolchain.from_env())

# `+` concatenates the two designs into one combined list of simulations
# (1200 training + a handful of held-out ones). prepare() then creates one
# folder per simulation under runs/paper/, containing everything penEasy
# needs: penEasy.in, the .geo geometry, the .spc spectrum, and the mat.mat /
# Aire.mat material files (built for you by calling material.exe).
study.prepare(training + held_out)
print(study.status())   # e.g. Counter({'pending': 1220})

# Run penEasy.x for every pending folder, up to 10 processes at a time. This
# is resumable: if it gets interrupted, running the script again picks up
# where it left off instead of starting over or re-running finished work.
study.run(max_concurrent=10)

# Harvest every finished simulation into one CSV row per simulation
# (composition, thickness, beam, and the computed transmission ratio).
study.collect(output="runs/paper/results.csv")
