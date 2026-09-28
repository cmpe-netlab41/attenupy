"""Gd/Pb thickness sweep at 100 kV: 20 simulations.

A minimal, complete example: two single-element shielding compositions (Gd
and Pb, each at 70 wt% in a C2H4 matrix) swept over ten thicknesses at one
beam. If you're new to attenupy, start here -- it's the shortest script that
still goes through all four stages (design, prepare, run, collect). This
reproduces (and lets you regenerate) the historical `simulations_copy/100kV`
folders:

    RQR8_100_1_03_oneL015_Gd_0.70 ... RQR8_100_1_03_oneL150_Pb_0.70

Run it with your own penEasy.x and material.exe:

    export ATTENUPY_PENEASY=/path/to/penEasy.x
    export ATTENUPY_MATERIAL=/path/to/pendbase/material.exe
    python examples/03_gd_pb_sweep.py

The same design, without writing any Python, is in
examples/gd_pb_sweep.toml -- see the comments in that file for the
equivalent `attenupy prepare` command.
"""

import attenupy as atp

# Print progress ("Prepared 20 simulation folder(s)...", "Launched 20...",
# etc.) to the terminal as the script runs.
atp.log_to_console()

# ThicknessSweep pairs every mixture in the list with every thickness from
# thickness_range(...), for every beam in beams=[...]: here that's
# 2 mixtures x 10 thicknesses x 1 beam = 20 simulations in total.
design = atp.ThicknessSweep(
    [atp.Mixture({"Gd": 0.70}, polymer_fraction=0.30),   # 70% Gd, 30% C2H4
     atp.Mixture({"Pb": 0.70}, polymer_fraction=0.30)],  # 70% Pb, 30% C2H4
    atp.thickness_range(0.15, 1.5, 0.15),   # 0.15, 0.30, ..., 1.50 mm
    beams=["RQR8_100"],                    # 100 kV
)

# Toolchain.from_env() reads penEasy.x / material.exe from the
# ATTENUPY_PENEASY / ATTENUPY_MATERIAL environment variables set above.
study = atp.Study("runs/gd_pb", toolchain=atp.Toolchain.from_env())

# Create the 20 folders: each gets penEasy.in, the .geo slab geometry, the
# .spc beam spectrum, and mat.mat / Aire.mat (built by calling material.exe).
study.prepare(design)
print(study.status())   # Counter({'pending': 20})

# Run penEasy.x for all 20 folders, up to 10 processes at once. Running the
# script again would resume instead of re-running finished simulations.
study.run(max_concurrent=10)

# Harvest the finished simulations into one CSV row per simulation.
study.collect(output="runs/gd_pb/results.csv")
