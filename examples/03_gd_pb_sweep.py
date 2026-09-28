"""Gd/Pb thickness sweep at 100 kV: 20 simulations.

Gd and Pb at 70 wt% in C2H4 (30 wt% polymer), RQR8_100 beam, 10 thicknesses
from 0.15 to 1.5 mm:

    RQR8_100_1_03_oneL015_Gd_0.70 ... RQR8_100_1_03_oneL150_Pb_0.70

    export ATTENUPY_PENEASY=/path/to/penEasy.x
    export ATTENUPY_MATERIAL=/path/to/pendbase/material.exe
    python examples/03_gd_pb_sweep.py

Same design as a file for the CLI: examples/gd_pb_sweep.toml.
"""

import attenupy as atp

atp.log_to_console()

design = atp.ThicknessSweep(
    [atp.Mixture({"Gd": 0.70}, polymer_fraction=0.30),
     atp.Mixture({"Pb": 0.70}, polymer_fraction=0.30)],
    atp.thickness_range(0.15, 1.5, 0.15),   # 0.15, 0.30, ..., 1.50 mm
    beams=["RQR8_100"],                    # 100 kV
)

study = atp.Study("runs/gd_pb", toolchain=atp.Toolchain.from_env())
study.prepare(design)   # 20 folders: penEasy.in, .geo, .spc, mat.mat, Aire.mat
print(study.status())

study.run(max_concurrent=10)
study.collect(output="runs/gd_pb/results.csv")
