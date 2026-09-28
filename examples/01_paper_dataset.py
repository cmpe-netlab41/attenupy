"""Build the paper's dataset design end to end.

    export ATTENUPY_PENEASY=/path/to/penEasy.x
    export ATTENUPY_MATERIAL=/path/to/pendbase/material.exe
    python examples/01_paper_dataset.py
"""

import attenupy as atp

atp.log_to_console()

SHIELDERS = ["Bi2O3", "BaSO4", "BN", "B4C", "WO3", "Sb", "Sn", "Bi", "W", "Cu", "B", "Ba"]
THICKNESSES = atp.thickness_range(0.15, 1.5, 0.15)  # mm

# Training set: 300 random composites per beam; 1-3 shielders in a C2H4 matrix,
# polymer fraction ~ U{0.2, 0.3, 0.4}, the rest split by a flat Dirichlet;
# each composition appears once, at one random thickness.
training = atp.RandomSampler(
    n=300,
    beams=["N30_030", "RQR6_080", "RQR8_100", "x120kV"],
    constituents=SHIELDERS,
    n_constituents=[1, 2, 3],
    polymer_fractions=[0.20, 0.30, 0.40],
    thicknesses_mm=THICKNESSES,
    unique="mixture",
    seed=2026,
)

# Held-out curves: fixed compositions at every thickness.
held_out = atp.ThicknessSweep(
    [atp.Mixture({"Bi2O3": 0.40, "BaSO4": 0.20}, polymer_fraction=0.40),
     atp.Mixture({"Bi": 0.35, "W": 0.35}, polymer_fraction=0.30)],
    THICKNESSES,
    beams=["RQR6_080"],
)

study = atp.Study("runs/paper", toolchain=atp.Toolchain.from_env())
study.prepare(training + held_out)   # folders, inputs, geometry, spectra, materials
print(study.status())

study.run(max_concurrent=10)
study.collect(output="runs/paper/results.csv")
