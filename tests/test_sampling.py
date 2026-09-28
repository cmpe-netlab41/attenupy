import os
import random
import unittest
from collections import Counter

import attenupy as atp
from helpers import TempDirTestCase

POOL = ["Bi2O3", "BaSO4", "BN", "B4C", "WO3", "Sb", "Sn", "Bi", "W", "Cu", "B", "Ba"]
THICK = atp.thickness_range(0.15, 1.5, 0.15)


def total(m):
    return m.polymer_fraction + sum(f for _, f in m.fractions)


class RandomSamplerTest(unittest.TestCase):
    def test_counts_sizes_and_sums(self):
        specs = atp.RandomSampler(n=120, beams=["RQR6_080", "RQR8_100"], constituents=POOL,
                                 n_constituents=[1, 2, 3], thicknesses_mm=THICK,
                                 seed=1).sample()
        self.assertEqual(Counter(s.beam for s in specs), {"RQR6_080": 120, "RQR8_100": 120})
        self.assertEqual(Counter(s.mixture.n_constituents for s in specs), {1: 80, 2: 80, 3: 80})
        for s in specs:
            self.assertAlmostEqual(total(s.mixture), 1.0, places=12)
            self.assertIn(s.mixture.polymer_fraction, (0.2, 0.3, 0.4))
            self.assertIn(s.thickness_mm, THICK)
            self.assertTrue(all(f >= 0.01 for _, f in s.mixture.fractions))

    def test_reproducible_and_seed_sensitive(self):
        make = lambda seed: atp.RandomSampler(n=20, beams=["80"], constituents=POOL,
                                             thicknesses_mm=THICK, seed=seed).sample()
        self.assertEqual(make(3), make(3))
        self.assertNotEqual(make(3), make(4))

    def test_size_weights(self):
        specs = atp.RandomSampler(n=100, beams=["80"], constituents=POOL,
                                 n_constituents={2: 0.7, 3: 0.3}, thicknesses_mm=THICK,
                                 seed=0).sample()
        self.assertEqual(Counter(s.mixture.n_constituents for s in specs), {2: 70, 3: 30})

    def test_balanced_combinations(self):
        # 66 pairs from 12; 132 draws -> every pair exactly twice
        specs = atp.RandomSampler(n=132, beams=["80"], constituents=POOL, n_constituents=2,
                                 thicknesses_mm=THICK, seed=0).sample()
        self.assertEqual(set(Counter(tuple(s.mixture.names) for s in specs).values()), {2})

    def test_require_any(self):
        specs = atp.RandomSampler(n=50, beams=["80"], constituents=POOL, n_constituents=2,
                                 require_any=["B", "Ba"], thicknesses_mm=THICK, seed=0).sample()
        self.assertTrue(all({"B", "Ba"} & set(s.mixture.names) for s in specs))

    def test_uniform_mode_and_full_precision(self):
        specs = atp.RandomSampler(n=30, beams=["80"], constituents=POOL, n_constituents=[2, 3],
                                 combinations="uniform", decimals=None,
                                 polymer_fractions=[0.25], thicknesses_mm=[0.5],
                                 seed=0).sample()
        self.assertEqual(len(specs), 30)
        self.assertTrue(any(len(repr(f)) > 5 for s in specs for _, f in s.mixture.fractions))

    def test_counts_per_beam_dict(self):
        specs = atp.RandomSampler(n={"80": 3, "120": 5}, constituents=POOL,
                                 thicknesses_mm=THICK, seed=0).sample()
        self.assertEqual(Counter(s.beam for s in specs), {"80": 3, "120": 5})

    def test_unique_specs_by_default(self):
        # 1-constituent space per beam: 12 x 3 polymer levels x 10 thicknesses = 360
        specs = atp.RandomSampler(n=200, beams=["80"], constituents=POOL, n_constituents=1,
                                 thicknesses_mm=THICK, seed=0).sample()
        self.assertEqual(len(set(specs)), 200)

    def test_unique_mixtures(self):
        specs = atp.RandomSampler(n=40, beams=["80", "100"], constituents=POOL,
                                 n_constituents=2, thicknesses_mm=THICK, unique="mixture",
                                 seed=0).sample()
        self.assertEqual(len({s.mixture for s in specs}), 80)

    def test_warns_when_space_too_small(self):
        with self.assertLogs("attenupy.sampling", "WARNING"):
            atp.RandomSampler(n=10, beams=["80"], constituents=["W"], polymer_fractions=[0.3],
                             n_constituents=1, thicknesses_mm=[1.0], max_redraws=3).sample()

    def test_bad_configs(self):
        with self.assertRaises(ValueError):
            atp.RandomMixtures(["W"], n_constituents=2)
        with self.assertRaises(ValueError):
            atp.RandomMixtures(POOL, require_any=["Pb"])
        with self.assertRaises(ValueError):
            atp.RandomMixtures(POOL, polymer_fractions=[1.2])

    def test_custom_mixture_source(self):
        class OnlyTungsten:
            def draw(self, count, rng):
                return [atp.Mixture({"W": 0.6}, 0.4)] * count

        specs = atp.RandomSampler(n=4, beams=["80"], mixtures=OnlyTungsten(),
                                 thicknesses_mm=[1.0]).sample()
        self.assertEqual({s.mixture.label for s in specs}, {"W_0.60"})


class SplitterTest(unittest.TestCase):
    def test_dirichlet_bias_and_floor(self):
        rng = random.Random(0)
        d = atp.Dirichlet({"B": 2.5}, min_share=0.1)
        shares = [d.split(["B", "W"], 0.7, rng) for _ in range(4000)]
        mean_b = sum(s[0] for s in shares) / len(shares) / 0.7
        self.assertAlmostEqual(mean_b, 0.1 + 0.8 * 2.5 / 3.5, delta=0.02)
        self.assertTrue(all(min(s) >= 0.07 - 1e-12 for s in shares))
        self.assertTrue(all(abs(sum(s) - 0.7) < 1e-12 for s in shares))

    def test_round_preserving_sum(self):
        rng = random.Random(0)
        for _ in range(500):
            k = rng.randint(1, 5)
            vals = [rng.random() for _ in range(k)]
            scale = 0.6 / sum(vals)
            out = atp.round_preserving_sum([v * scale for v in vals], 2)
            self.assertAlmostEqual(sum(out), 0.6, places=9)
            self.assertTrue(all(x >= 0.01 for x in out))
        self.assertEqual(atp.round_preserving_sum([1 / 3] * 3, 2), [0.34, 0.33, 0.33])


class GridAndSweepTest(TempDirTestCase):
    def test_sweep(self):
        mix = atp.Mixture({"Bi2O3": 0.4, "BaSO4": 0.2}, 0.4)
        specs = atp.ThicknessSweep(mix, THICK, beams=["80", "100"]).sample()
        self.assertEqual(len(specs), 20)
        self.assertEqual(sorted({s.thickness_mm for s in specs}), THICK)

    def test_sweep_random_mixtures_and_chain(self):
        mixes = atp.RandomMixtures(POOL, n_constituents=2).sample(6, seed=1)
        sweep = atp.ThicknessSweep(mixes, THICK, beams=["80"])
        combined = sweep + atp.GridSampler(["120"], mixes[:1], [0.5])
        self.assertEqual(len(combined.sample()), 61)

    def test_thickness_range(self):
        self.assertEqual(atp.thickness_range(0.15, 1.5, 0.15),
                         [0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.05, 1.2, 1.35, 1.5])

    def test_wide_csv(self):
        path = os.path.join(self.tmp, "rows.csv")
        header = ([f"molecule_{i}" for i in range(1, 11)] + [f"ratio_{i}" for i in range(1, 11)]
                  + ["polymer", "polymer_ratio", "geo_number_mm", "spc"])
        row = (["Sb", "Cu"] + ["-1"] * 8 + ["20.0", "60.0"] + ["-1.0"] * 8
               + ["C2H4", "0.2", "15", "120"])
        with open(path, "w") as f:
            f.write(",".join(header) + "\n" + ",".join(row) + "\n")
        (s,) = atp.CSVSampler(path).sample()
        self.assertEqual((s.beam, s.thickness_mm), ("120", 0.15))
        self.assertEqual(s.mixture.as_dict(), {"Cu": 0.6, "Sb": 0.2})

    def test_csv_custom_parser_multiple_specs_per_row(self):
        path = os.path.join(self.tmp, "70kv.csv")
        with open(path, "w") as f:
            f.write("PU,Bi2O3,t1,t2\n40,60,0.5,1.0\n")

        def parse(row):
            mix = atp.Mixture({"Bi2O3": float(row["Bi2O3"]) / 100}, float(row["PU"]) / 100)
            return [atp.SimulationSpec("RQR5_070", mix, float(row[c])) for c in ("t1", "t2")]

        self.assertEqual([s.thickness_mm for s in atp.CSVSampler(path, parse)], [0.5, 1.0])


if __name__ == "__main__":
    unittest.main()
