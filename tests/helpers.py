"""Test doubles: a stand-in penEasy.x and a material builder without material.exe."""

import os
import stat
import sys
import tempfile
import unittest

DATA = os.path.join(os.path.dirname(__file__), "data")

# Reads penEasy.in from stdin like penEasy, writes an energy-deposition tally
# where the detector (material 2) sees I0 * exp(-mu * t), and signs off.
# FAKE_PENEASY_FAIL=1 makes it crash instead.
FAKE_PENEASY = r'''#!{python}
import math, os, re, sys
deck = sys.stdin.read()
if os.environ.get("FAKE_PENEASY_FAIL"):
    print("simulated crash"); sys.exit(3)
geo = re.search(r"^ (\S+\.geo) ", deck, re.M).group(1)
slab = re.search(r"^  1   (\S+)", deck, re.M).group(1)
t_mm = float(open(geo).readline().split()[0])
i0 = 0.8
det = i0 if slab == "Aire.mat" else i0 * math.exp(-2.0 * t_mm)
with open("tallyEnergyDeposition.dat", "w") as f:
    f.write("# Material : Energy (eV/hist) : +-2sigma\n")
    f.write(f"   1  1.0E+03  1.0E-01\n   2  {{det:.5E}}  5.0E-03\n   3  1.0E+03  1.0E-01\n")
print("Program ended")
print(" ")
print("Have a nice day.")
'''


def make_fake_peneasy(directory):
    path = os.path.join(directory, "penEasy.x")
    with open(path, "w") as f:
        f.write(FAKE_PENEASY.format(python=sys.executable))
    os.chmod(path, os.stat(path).st_mode | stat.S_IEXEC)
    return path


class FakeMaterials:
    """Writes .mat stubs carrying the mixture density, and counts calls."""

    def __init__(self):
        self.built = []
        self.air_built = 0

    def build(self, mixture, dest, registry=None):
        self.built.append(mixture)
        with open(dest, "w") as f:
            f.write(f" Material: {mixture.label}\n")
            f.write(f" Mass density = {mixture.density(registry):.8E} g/cm**3\n")
            f.write(" Mean excitation energy = 1.23450000E+02 eV\n")

    def build_air(self, dest):
        self.air_built += 1
        with open(dest, "w") as f:
            f.write(" Material: AIR\n Mass density = 1.20479000E-03 g/cm**3\n")


class TempDirTestCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()
