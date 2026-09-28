"""Generated input files match known-good reference files exactly."""

import os
import unittest

import attenupy as atp
from helpers import DATA, TempDirTestCase


def spec(mm, beam="RQR8_100"):
    return atp.SimulationSpec(beam, atp.Mixture({"Gd": 0.70}, 0.30), mm)


class GeometryTest(unittest.TestCase):
    def test_matches_reference_files(self):
        for name, mm in [("oneL015", 0.15), ("oneL250", 2.5)]:
            with open(os.path.join(DATA, f"{name}.geo")) as f:
                self.assertEqual(atp.SlabGeometry().render(spec(mm)), f.read(), name)

    def test_matches_hand_edited_file_up_to_trailing_spaces(self):
        with open(os.path.join(DATA, "oneL030.geo")) as f:
            ref = [line.rstrip() for line in f]
        out = [line.rstrip() for line in atp.SlabGeometry().render(spec(0.30)).splitlines()]
        self.assertEqual(out, ref)

    def test_filename_and_overlap_check(self):
        g = atp.SlabGeometry()
        self.assertEqual(g.filename(spec(0.15)), "oneL015.geo")
        self.assertEqual(g.filename(spec(0.125)), "slab_0.125mm.geo")
        with self.assertRaises(ValueError):
            g.validate(spec(10.0))  # would reach the detector at 1 cm
        atp.SlabGeometry(detector_z_cm=102).validate(spec(10.0))

    def test_custom_dimensions(self):
        text = atp.SlabGeometry(slab_half_width_cm=5, detector_z_cm=110).render(spec(1.0))
        self.assertIn("10x10 cm2", text)
        self.assertIn("A0=( 0.500000000000000E+01,   0)", text)
        self.assertIn("Z-SHIFT=( 1.100000000000000E+02,   0)", text)
        self.assertIn("at 110 cm from the source, 10 cm from", text)

    def test_beam_must_fit_on_the_slab(self):
        atp.SlabGeometry().check_beam(2.57)  # default setup: 4.49 cm cone on a 4.5 cm slab
        with self.assertRaisesRegex(ValueError, "unshielded"):
            atp.SlabGeometry(slab_half_width_cm=3.0).check_beam(2.57)
        with self.assertRaisesRegex(ValueError, "unshielded"):
            atp.SlabGeometry().check_beam(5.0)

    def test_seeds_are_validated(self):
        mix = atp.Mixture({"Gd": 0.70}, 0.30)
        self.assertEqual(atp.SimulationSpec("RQR8_100", mix, 1.0, seed=7).seed, (7, 7))
        for bad in [(0, 0), (-1, -1), 0, (1, 2147483399)]:
            with self.assertRaisesRegex(ValueError, "seeds", msg=bad):
                atp.SimulationSpec("RQR8_100", mix, 1.0, seed=bad)


class DeckTest(TempDirTestCase):
    def test_default_deck_matches_reference_file(self):
        with open(os.path.join(DATA, "penEasy_RQR8_100_oneL015.in")) as f:
            ref = f.read()
        self.assertEqual(atp.PeneasyDeck().render("RQR8_100.spc", "oneL015.geo"), ref)

    def test_knobs_and_reference_material(self):
        deck = atp.PeneasyDeck(num_histories="1.0e6").render(
            "A.spc", "g.geo", seed=(7, 9), slab_material="Aire.mat")
        self.assertIn(" 1.0e6             NUMBER OF HISTORIES", deck)
        self.assertIn(" 7  9              INITIAL RANDOM SEEDS", deck)
        self.assertIn("  1   Aire.mat          100.0e3", deck)
        with self.assertRaises(ValueError):
            atp.PeneasyDeck().render("A.spc", "g.geo", slab_material="x" * 21)

    def test_custom_template_file(self):
        path = os.path.join(self.tmp, "t.in")
        with open(path, "w") as f:
            f.write("{spectrum_file} {geometry_file} {thickness_mm} {{literal}}\n")
        deck = atp.PeneasyDeck.from_template_file(path)
        self.assertEqual(deck.render("a.spc", "b.geo", thickness_mm=0.5), "a.spc b.geo 0.5 {literal}\n")


class SpectrumTest(TempDirTestCase):
    def test_bundled_spectra(self):
        self.assertEqual(atp.bundled_spectra(),
                         ["N30_030", "RQR5_070", "RQR6_080", "RQR8_100", "x120kV"])
        e, w = atp.bundled_spectrum("RQR6_080").read()
        self.assertEqual(e[-1], 80000.0)
        self.assertTrue(all(b > a for a, b in zip(e, e[1:])))

    def test_tabulated_roundtrip(self):
        path = os.path.join(self.tmp, "x.spc")
        atp.TabulatedSpectrum([10e3, 20e3, 30e3], [1, 2, 0.5], header="toy").write(path)
        self.assertEqual(atp.read_spectrum(path), ([10e3, 20e3, 30e3], [1.0, 2.0, 0.5]))
        with open(path) as f:  # last bin as wide as the previous one, not zero-width
            self.assertEqual(f.read().splitlines()[-1], "40000\t-1")

    def test_from_points_uses_trapezoids(self):
        # penEasy manual, Appendix B: P_k = (f_k + f_k+1)/2 * (E_k+1 - E_k)
        spc = atp.TabulatedSpectrum.from_points([10e3, 20e3, 40e3], [2.0, 4.0, 0.0])
        self.assertEqual(spc.energies_eV, (10e3, 20e3))
        self.assertEqual(spc.weights, (3.0e4, 4.0e4))
        self.assertEqual(spc.end_energy_eV, 40e3)

    def test_monoenergetic_line(self):
        path = os.path.join(self.tmp, "line.spc")
        atp.write_bigs_spectrum(path, [59.5e3], [1.0])
        self.assertEqual(atp.read_spectrum(path), ([59.5e3], [1.0]))
        self.assertEqual(atp.check_bigs_spectrum(path), 1)

    def test_equal_sources_compare_equal(self):
        self.assertEqual(atp.bundled_spectrum("RQR6_080"), atp.bundled_spectrum("RQR6_080"))


class ChemistryTest(unittest.TestCase):
    def test_parse_formula(self):
        self.assertEqual(atp.parse_formula("Bi2O3"), {"Bi": 2, "O": 3})
        self.assertEqual(atp.parse_formula("Mg(OH)2"), {"Mg": 1, "O": 2, "H": 2})
        self.assertEqual(atp.parse_formula("Ca10(PO4)6(OH)2"), {"Ca": 10, "P": 6, "O": 26, "H": 2})
        for bad in ["", "Xx2", "bi2", "Bi2O3)", "(Bi2", "2Bi"]:
            with self.assertRaises(ValueError, msg=bad):
                atp.parse_formula(bad)

    def test_density_matches_reference_mat_file(self):
        # mat.mat of RQR8_100_1_03_oneL015_Gd_0.70 reports 2.47281700 g/cm**3
        self.assertAlmostEqual(atp.Mixture({"Gd": 0.70}, 0.30).density(), 2.472817, places=6)

    def test_registry(self):
        reg = atp.ConstituentRegistry()
        reg.register("PU", density=1.2, formula="C3H8N2O")
        self.assertEqual(reg["PU"].formula, "C3H8N2O")
        with self.assertRaises(ValueError):
            reg.register("PU", density=1.3)
        with self.assertRaisesRegex(KeyError, "register_constituent"):
            reg["Unobtainium"]
        with self.assertRaises(ValueError):
            reg.register("my_thing", density=1.0, formula="C")  # '_' is a separator


class BeamTest(unittest.TestCase):
    def test_aliases(self):
        for alias in ["RQR6_080", "80", "80kV"]:
            self.assertEqual(atp.DEFAULT_BEAMS[alias].key, "RQR6_080")

    def test_ambiguous_alias(self):
        reg = atp.DEFAULT_BEAMS.copy()
        reg.register("W80", 80, atp.bundled_spectrum("RQR6_080"), folder="80kV_W")
        with self.assertRaisesRegex(KeyError, "ambiguous"):
            reg["80"]
        self.assertEqual(reg["80kV"].key, "RQR6_080")

    def test_reregistering_same_beam_is_fine(self):
        reg = atp.BeamRegistry()
        reg.register("RQR6_080", 80)
        reg.register("RQR6_080", 80)


class NamingTest(unittest.TestCase):
    def test_compact_roundtrip(self):
        name = "RQR6_080_2_04_oneL080_BaSO4_Bi2O3_0.20_0.40"
        s = atp.CompactNamer.parse(name)
        self.assertEqual(s.mixture.as_dict(), {"BaSO4": 0.2, "Bi2O3": 0.4})
        self.assertEqual(s.thickness_mm, 0.8)
        self.assertEqual(atp.CompactNamer().name(s), name)

    def test_compact_refuses_to_mislabel(self):
        s = atp.SimulationSpec("RQR6_080", atp.Mixture({"W": 0.75}, 0.25), 1.0)
        with self.assertRaisesRegex(ValueError, "DescriptiveNamer"):
            atp.CompactNamer().name(s)
        self.assertEqual(atp.DescriptiveNamer().name(s), "RQR6_080_n1_p25_t1mm_W-0.7500")


if __name__ == "__main__":
    unittest.main()
