"""Regression checks for rejecting incompatible/stale comparison tables."""
import importlib.util
import tempfile
import unittest
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location("comparison", ROOT / "scripts/figures/plot_model_comparison_figures.py")
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)

class ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.path = Path(self.temp.name) / "results.csv"
        self.rows = pd.DataFrame([
            dict(model=model, dic=i+100, waic=i+200, mae=2.0, rmse=3.0, wape=.5, r2=.2,
                 evaluation_protocol="common_observations_2023_january_v1",
                 evaluation_signature="shared", test_n=100, test_start="2023-01-02",
                 test_end="2023-01-23", fit_signature="one" if i < 5 else "two",
                 fit_n=2000, criteria_group="fit_1" if i < 5 else "fit_2", run_id="run")
            for i, model in enumerate(comparison.MODEL_ORDER)])
    def tearDown(self):
        self.temp.cleanup()
    def load(self):
        self.rows.to_csv(self.path, index=False)
        return comparison.load_model_results(self.path)
    def test_grouped_criteria(self):
        result = self.load()
        self.assertEqual(result.loc[5, "delta_waic"], 0)
        self.assertEqual(result.loc[5, "waic_rank_recomputed"], 1)
    def test_different_sample(self):
        self.rows.loc[1, "evaluation_signature"] = "different"
        with self.assertRaisesRegex(ValueError, "same evaluation contract"):
            self.load()
    def test_missing_model(self):
        self.rows = self.rows.iloc[:-1]
        with self.assertRaisesRegex(ValueError, "17 models"):
            self.load()
    def test_invalid_metric(self):
        self.rows.loc[1, "rmse"] = float("inf")
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            self.load()
    def test_legacy_table(self):
        self.rows = self.rows.drop(columns="evaluation_signature")
        with self.assertRaisesRegex(ValueError, "Missing required"):
            self.load()

if __name__ == "__main__":
    unittest.main()
