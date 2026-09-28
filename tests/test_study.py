import csv
import json
import logging
import math
import os
import shutil
import unittest
from unittest import mock

import attenupy as atp
from attenupy.cli import main as cli_main
from helpers import FakeMaterials, TempDirTestCase, make_fake_peneasy

MIX = atp.Mixture({"Bi2O3": 0.40, "BaSO4": 0.20}, 0.40)
REAL_MATERIAL = os.environ.get("ATTENUPY_MATERIAL") or os.path.join(
    os.path.dirname(__file__), "..", "..", "new_simulations_can", "penelope_peneasy",
    "penelope", "pendbase", "material.exe")


class StudyTest(TempDirTestCase):
    def make_study(self, **kwargs):
        self.fake = FakeMaterials()
        self.exe = make_fake_peneasy(self.tmp)
        return atp.Study(os.path.join(self.tmp, "study"), materials=self.fake,
                        toolchain=atp.Toolchain(peneasy=self.exe), **kwargs)

    def test_prepare_writes_everything(self):
        study = self.make_study()
        (path,) = study.prepare([atp.SimulationSpec("80", MIX, 0.8)])
        self.assertEqual(os.path.relpath(path, study.root),
                         os.path.join("80kV", "RQR6_080_2_04_oneL080_BaSO4_Bi2O3_0.20_0.40"))
        self.assertEqual(sorted(f for f in os.listdir(path) if not f.startswith(".")),
                         ["Aire.mat", "RQR6_080.spc", "mat.mat", "oneL080.geo",
                          "penEasy.in", "spec.json"])
        saved = atp.SimulationSpec.load(path)
        self.assertEqual(saved, atp.SimulationSpec("RQR6_080", MIX, 0.8))
        with open(os.path.join(path, "penEasy.in")) as f:
            deck = f.read()
        self.assertIn(" RQR6_080.spc ", deck)
        self.assertIn(" oneL080.geo ", deck)
        self.assertTrue(os.path.isfile(os.path.join(study.state_dir, "study.json")))

    def test_inputs_only_then_materials(self):
        study = self.make_study()
        study.prepare([atp.SimulationSpec("80", MIX, 0.8)], build_materials=False)
        self.assertEqual(study.status(), {"unprepared": 1})
        self.assertEqual(study.run(poll_seconds=0.05).launched, 0)
        study.build_materials()
        self.assertEqual(study.status(), {"pending": 1})

    def test_empty_material_counts_as_missing(self):
        # penEasy creates an empty file when a material is missing
        study = self.make_study()
        (path,) = study.prepare([atp.SimulationSpec("80", MIX, 0.8)])
        open(os.path.join(path, "mat.mat"), "w").close()
        self.assertEqual(study.status(), {"unprepared": 1})
        study.build_materials()
        self.assertEqual(study.status(), {"pending": 1})

    def test_plan_rejects_a_beam_wider_than_the_slab(self):
        study = self.make_study(geometry=atp.SlabGeometry(slab_half_width_cm=3.0))
        with self.assertRaisesRegex(ValueError, "unshielded"):
            study.prepare([atp.SimulationSpec("80", MIX, 0.8)])
        self.assertFalse(os.path.exists(study.root))

    def test_sweep_builds_each_material_once(self):
        study = self.make_study()
        study.prepare(atp.ThicknessSweep(MIX, atp.thickness_range(0.15, 1.5, 0.15), beams=["80"]))
        self.assertEqual(len(self.fake.built), 1)
        self.assertEqual(self.fake.air_built, 1)
        study.prepare(atp.ThicknessSweep(MIX, [0.15, 1.65], beams=["100"]))
        self.assertEqual(len(self.fake.built), 1)  # cache hit across calls

    def test_invalid_design_writes_nothing(self):
        study = self.make_study()
        bad = [atp.SimulationSpec("80", MIX, 0.8),
               atp.SimulationSpec("999kV", MIX, 0.8),
               atp.SimulationSpec("80", atp.Mixture({"Unobtainium": 0.6}, 0.4), 0.8),
               atp.SimulationSpec("80", MIX, 12.0)]
        with self.assertRaises(ValueError) as cm:
            study.prepare(bad)
        msg = str(cm.exception)
        for fragment in ("999kV", "Unobtainium", "overlap"):
            self.assertIn(fragment, msg)
        self.assertFalse(os.path.exists(study.root))

    def test_duplicates_and_collisions(self):
        study = self.make_study()
        a = atp.SimulationSpec("80", atp.Mixture({"W": 0.6}, 0.4), 0.8)
        b = atp.SimulationSpec("80", atp.Mixture({"W": 0.6}, 0.4), 0.8, seed=(5, 5))
        paths = study.prepare([a, a, b])
        self.assertEqual(len(paths), 2)
        self.assertTrue(paths[1].endswith("__2"))
        # preparing again maps to the same folders
        self.assertEqual(study.prepare([a, b]), paths)

    def test_full_cycle_run_collect(self):
        study = self.make_study()
        study.prepare(atp.ThicknessSweep(MIX, [0.5, 1.0], beams=["80"]))
        self.assertEqual(study.status(), {"pending": 2})
        summary = study.run(max_concurrent=2, poll_seconds=0.05)
        self.assertEqual((summary.launched, summary.succeeded), (2, 2))
        self.assertEqual(study.status(), {"done": 2})
        self.assertEqual(study.run(poll_seconds=0.05).launched, 0)  # resumable

        out = os.path.join(self.tmp, "r.csv")
        rows = study.collect(output=out)
        self.assertEqual(len(rows), 2)
        i0 = atp.DEFAULT_BEAMS["80"].reference_energy
        for row in rows:
            t = row["geo_number_mm"] / 100
            self.assertAlmostEqual(row["energy_ratio"], 0.8 * math.exp(-2 * t) / i0, places=4)
            self.assertEqual((row["molecule_1"], row["ratio_1"], row["molecule_3"]),
                             ("BaSO4", 20.0, "-1"))
        with open(out) as f:
            header = next(csv.reader(f))
        self.assertEqual(header[-7:], ["polymer", "polymer_ratio", "geo_number_mm", "spc",
                                       "mass_den", "mean_excit", "energy_ratio"])

    def test_measured_reference_takes_precedence(self):
        study = self.make_study()
        study.prepare([atp.SimulationSpec("80", MIX, 1.0)])
        study.prepare_references()
        self.assertEqual(study.status(), {"pending": 2})
        study.run(poll_seconds=0.05)
        self.assertEqual(study.measured_references(), {"RQR6_080": 0.8})
        (row,) = study.collect()
        self.assertAlmostEqual(row["energy_ratio"], math.exp(-2.0), places=4)

    def test_new_beam_without_reference_energy(self):
        beams = atp.DEFAULT_BEAMS.copy()
        beams.register("W60", 60, atp.TabulatedSpectrum([10e3, 30e3, 60e3], [1, 3, 1]))
        study = self.make_study(beams=beams)
        (path,) = study.prepare([atp.SimulationSpec("W60", MIX, 0.5)])
        self.assertTrue(os.path.isfile(os.path.join(path, "W60.spc")))
        study.run(poll_seconds=0.05)
        (row,) = study.collect()
        self.assertIsNone(row["energy_ratio"])  # no I0 until a reference run exists

    def test_failures_are_recorded_and_retryable(self):
        study = self.make_study()
        study.prepare([atp.SimulationSpec("80", MIX, 0.5)])
        with mock.patch.dict(os.environ, {"FAKE_PENEASY_FAIL": "1"}):
            summary = study.run(poll_seconds=0.05)
        self.assertEqual(len(summary.failed), 1)
        self.assertEqual(study.status(), {"failed": 1})
        self.assertEqual(study.run(poll_seconds=0.05).launched, 0)
        self.assertEqual(study.run(retry_failed=True, poll_seconds=0.05).succeeded, 1)

    def test_claims_block_other_runners(self):
        study = self.make_study()
        (path,) = study.prepare([atp.SimulationSpec("80", MIX, 0.5)])
        os.makedirs(os.path.join(path, ".attenupy"), exist_ok=True)
        with open(os.path.join(path, ".attenupy", "claim"), "w") as f:
            json.dump({"host": "some-other-node", "pid": 1}, f)
        self.assertEqual(atp.sim_state(path), "running")
        self.assertEqual(study.run(poll_seconds=0.05).launched, 0)
        with open(os.path.join(path, ".attenupy", "claim"), "w") as f:
            json.dump({"host": __import__("socket").gethostname(), "pid": 2 ** 22 + 12345}, f)
        self.assertEqual(atp.sim_state(path), "pending")  # dead local process -> stale

    def test_no_rerun_when_another_runner_finished_it(self):
        study = self.make_study()
        (path,) = study.prepare([atp.SimulationSpec("80", MIX, 0.5)])
        study.run(poll_seconds=0.05)
        with open(os.path.join(path, "penEasy.out")) as f:
            first = f.read()
        # simulate a runner whose to-do list was selected before the run finished
        with mock.patch("attenupy.runner.select", return_value=[path]):
            summary = atp.LocalRunner(self.exe, poll_seconds=0.05).run([path])
        self.assertEqual((summary.launched, summary.skipped), (0, 1))
        self.assertFalse(os.path.exists(os.path.join(path, ".attenupy", "claim")))
        with open(os.path.join(path, "penEasy.out")) as f:
            self.assertEqual(f.read(), first)

    def test_stale_material_is_flagged(self):
        study = self.make_study()
        (path,) = study.prepare([atp.SimulationSpec("80", MIX, 0.5)])
        study.run(poll_seconds=0.05)
        self.fake.build(atp.Mixture({"BN": 0.8}, 0.2), os.path.join(path, "mat.mat"))
        ctx = atp.CollectContext()
        atp.collect(study.simulations(), ctx=ctx)
        self.assertEqual(len(ctx.density_mismatches), 1)

    def test_custom_extractor(self):
        study = self.make_study()
        study.prepare([atp.SimulationSpec("80", MIX, 0.5)])
        study.run(poll_seconds=0.05)

        def det_sigma(sim, ctx):
            return {"det_2sigma": atp.read_energy_deposition(sim.file("tallyEnergyDeposition.dat"))[2][1]}

        (row,) = study.collect([atp.SpecColumns(), det_sigma])
        self.assertEqual((row["w_Bi2O3"], row["thickness_mm"], row["det_2sigma"]), (0.4, 0.5, 5e-3))

    def test_tree_without_spec_json(self):
        study = self.make_study()
        (path,) = study.prepare([atp.SimulationSpec("80", MIX, 0.5)])
        os.remove(os.path.join(path, "spec.json"))
        (sim,) = study.simulations()
        self.assertEqual(sim.spec, atp.SimulationSpec("RQR6_080", MIX, 0.5))


class CLITest(TempDirTestCase):
    def tearDown(self):
        logger = logging.getLogger("attenupy")  # the CLI turns console logging on
        logger.handlers.clear()
        logger.propagate = True
        logger.setLevel(logging.NOTSET)
        super().tearDown()

    def run_cli(self, *args):
        with self.assertRaises(SystemExit) as cm:
            cli_main(list(args))
        return cm.exception.code

    def test_prepare_dry_run(self):
        design = os.path.join(self.tmp, "design.json")
        with open(design, "w") as f:
            json.dump({
                "constituents": {"Gd2O3": 7.41},
                "design": [
                    {"type": "random", "n": 10, "beams": ["80"], "seed": 0,
                     "constituents": ["Gd2O3", "W", "Bi"], "n_constituents": [1, 2],
                     "thicknesses_mm": {"start": 0.15, "stop": 0.6, "step": 0.15},
                     "splitter": {"type": "dirichlet", "alpha": {"Gd2O3": 2.0}}},
                    {"type": "sweep", "beams": ["100"], "thicknesses_mm": [0.5, 1.0],
                     "mixtures": [{"fractions": {"W": 0.6}, "polymer_fraction": 0.4}]},
                ]}, f)
        root = os.path.join(self.tmp, "s")
        self.assertEqual(self.run_cli("prepare", root, "--design", design, "--dry-run"), 0)
        self.assertFalse(os.path.exists(root))
        self.assertEqual(len(list(atp.load_design(design))), 12)

    def test_doctor_reports_missing(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            self.assertEqual(self.run_cli("doctor"), 1)


@unittest.skipUnless(os.path.isfile(REAL_MATERIAL) and not atp.Toolchain(
    material=REAL_MATERIAL).material_problems(), "real material.exe not available")
class RealMaterialTest(TempDirTestCase):
    def test_builds_correct_material_and_air(self):
        tc = atp.Toolchain(material=REAL_MATERIAL)
        study = atp.Study(os.path.join(self.tmp, "s"), toolchain=tc)
        mix3 = atp.Mixture({"Bi2O3": 0.2, "BaSO4": 0.2, "WO3": 0.2}, 0.4)
        paths = study.prepare([atp.SimulationSpec("100", atp.Mixture({"Gd": 0.7}, 0.3), 0.15),
                               atp.SimulationSpec("100", mix3, 0.15)])
        for p, mix in zip(paths, [atp.Mixture({"Gd": 0.7}, 0.3), mix3]):
            self.assertAlmostEqual(atp.mass_density(os.path.join(p, "mat.mat")),
                                   mix.density(), places=5)
            self.assertAlmostEqual(atp.mass_density(os.path.join(p, "Aire.mat")), 1.20479e-3)
        # long 3-constituent paths used to overflow material.exe's 80-char limit
        self.assertGreater(len(os.path.join(paths[1], "mat.mat")), 80)


if __name__ == "__main__":
    unittest.main()
