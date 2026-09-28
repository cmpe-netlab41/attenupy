"""Spectrum files must survive penEasy's reader; SpekPy spectra are well formed."""

import importlib.util
import json
import os
import sys
import unittest
from unittest import mock

import attenupy as atp
from helpers import TempDirTestCase

HAVE_SPEKPY = importlib.util.find_spec("spekpy") is not None


class CheckBigsTest(TempDirTestCase):
    def write(self, text):
        path = os.path.join(self.tmp, "x.spc")
        with open(path, "w") as f:
            f.write(text)
        return path

    def test_bundled_spectra_are_valid(self):
        for name in atp.bundled_spectra():
            self.assertGreater(atp.check_bigs_spectrum(atp.bundled_spectrum(name).path), 50)

    def test_bin_limit(self):
        n = atp.spectra.MAX_BINS + 1
        path = self.write("".join(f"{1000 + i} 1\n" for i in range(n)) + f"{1000 + n} -1\n")
        with self.assertRaisesRegex(ValueError, "bins"):
            atp.check_bigs_spectrum(path)

    def test_rejects_what_penEasy_rejects(self):
        cases = {
            "Energy(keV)\tFluence\n1000\t1\n2000\t-1\n": "column 1",   # header without '#'
            " # indented comment\n1000 1\n2000 -1\n": "column 1",       # '#' not in column 1
            "1000\t1\n2000\t2\n": "negative weight",                    # no terminator
            "2000 1\n1000 2\n3000 -1\n": "decreasing",
            "1000 0\n2000 -1\n": "no positive weight",
        }
        for text, fragment in cases.items():
            with self.assertRaisesRegex(ValueError, fragment, msg=text):
                atp.check_bigs_spectrum(self.write(text))

    def test_accepts_comments_blanks_and_tabs(self):
        self.assertEqual(atp.check_bigs_spectrum(self.write("# c\n\n1000\t1\n 2000\t2\n3000 -1\n")), 2)

    def test_file_spectrum_validation_catches_bad_files_before_prepare(self):
        path = self.write("Energy Fluence\n1000 1\n2000 -1\n")
        beams = atp.BeamRegistry()
        beams.register("BAD", 2, path)
        study = atp.Study(os.path.join(self.tmp, "s"), beams=beams, materials=object())
        with self.assertRaisesRegex(ValueError, "column 1"):
            study.prepare([atp.SimulationSpec("BAD", atp.Mixture({"W": 0.6}, 0.4), 1.0)])
        self.assertFalse(os.path.exists(study.root))

    def test_missing_spekpy_gives_install_hint(self):
        with mock.patch.dict(sys.modules, {"spekpy": None}):
            with self.assertRaisesRegex(ImportError, "pip install spekpy"):
                atp.SpekpySpectrum(70).validate()


@unittest.skipUnless(HAVE_SPEKPY, "spekpy not installed")
class SpekpyTest(TempDirTestCase):
    SPEC = dict(kvp=70, th=12, filters=[("Al", 2.8), ("Be", 3.0)])

    def test_rqr5_file(self):
        path = os.path.join(self.tmp, "RQR5.spc")
        spc = atp.SpekpySpectrum(**self.SPEC)
        spc.write(path)
        self.assertEqual(atp.check_bigs_spectrum(path), 69)
        with open(path) as f:
            lines = f.read().splitlines()
        self.assertEqual(lines[-1], "70000\t-1")
        self.assertIn("HVL1 = 2.58", "\n".join(lines[:5]))
        self.assertAlmostEqual(spc.hvl()[0], 2.58, delta=0.01)  # IEC 61267 RQR5

        energies, weights = atp.read_spectrum(path)
        self.assertEqual(energies, [1000.0 * i for i in range(1, 70)])
        k, f = spc.spek().get_spectrum()
        for w, ref in zip(weights, f):
            self.assertAlmostEqual(w, ref, delta=abs(ref) * 1e-7)

    def test_half_kev_bins(self):
        edges, fluence, end = atp.SpekpySpectrum(30, th=18, dk=0.5, filters=[("Al", 4.0)]).table()
        self.assertEqual((edges[1] - edges[0], end), (500.0, 30000.0))

    def test_equal_specs_compare_equal(self):
        self.assertEqual(atp.SpekpySpectrum(**self.SPEC), atp.SpekpySpectrum(**self.SPEC))

    def test_design_file_beam(self):
        path = os.path.join(self.tmp, "d.json")
        with open(path, "w") as f:
            json.dump({
                "beams": [{"key": "RQR5_spekpy", "kvp": 70, "folder": "70kV_spekpy",
                           "spectrum": {"type": "spekpy", "th": 12,
                                        "filters": [["Al", 2.8], ["Be", 3.0]]}}],
                "design": [{"type": "sweep", "beams": ["RQR5_spekpy"], "thicknesses_mm": [0.5],
                            "mixtures": [{"fractions": {"W": 0.6}, "polymer_fraction": 0.4}]}],
            }, f)
        beams = atp.DEFAULT_BEAMS.copy()
        (spec,) = atp.load_design(path, beams=beams)
        self.assertEqual(beams[spec.beam].spectrum, atp.SpekpySpectrum(**self.SPEC))


if __name__ == "__main__":
    unittest.main()
