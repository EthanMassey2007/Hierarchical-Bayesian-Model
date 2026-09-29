#!/usr/bin/env python3
"""Write a concise manuscript-facing audit from a verified comparison run."""
import csv
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "outputs"

def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))

def main():
    status = read_csv(OUTPUT / "evaluation_status.csv")
    if len(status) != 1 or status[0]["status"] != "complete":
        raise SystemExit("Evaluation must complete successfully before reporting.")
    rows = read_csv(OUTPUT / "all_model_results_table.csv")
    if len(rows) != 17 or {r["run_id"] for r in rows} != {status[0]["run_id"]}:
        raise SystemExit("Missing models or stale aggregate table.")
    for field in ("test_n", "test_start", "test_end", "evaluation_signature"):
        if len({r[field] for r in rows}) != 1:
            raise SystemExit(f"Inconsistent common-sample field: {field}")
    run_dir = OUTPUT / "evaluation_runs" / rows[0]["run_id"]
    audit = read_csv(run_dir / "common_sample_audit.csv")
    cohort = read_csv(run_dir / "common_test_observations.csv")
    municipalities = len({r["municipio"] for r in cohort})
    lines = ["# Corrected common-sample evaluation", "",
             f"Run: `{rows[0]['run_id']}`", "",
             f"All 17 models were trained on available 2017–2022 observations and scored on "
             f"the same **{rows[0]['test_n']} municipality-week observations** "
             f"({municipalities} municipalities), dated **{rows[0]['test_start']} to {rows[0]['test_end']}**.", "",
             "The previous cross-model predictive metrics used unequal observation sets and must be replaced. "
             "Model-specific training samples still differ; this comparison evaluates fitted pipelines, "
             "rather than isolating the effect of a feature on an identical training cohort.", "",
             "## Replacement predictive metrics", "",
             "| Model | MAE | RMSE | WAPE | R² | N |", "|---|---:|---:|---:|---:|---:|"]
    for row in rows:
        lines.append(f"| {row['model'].replace('R_M','M')} | " + " | ".join(
            f"{float(row[k]):.4f}" for k in ("mae", "rmse", "wape", "r2")) + f" | {row['test_n']} |")
    best_wape = min(rows, key=lambda r: float(r["wape"]))
    best_rmse = min(rows, key=lambda r: float(r["rmse"]))
    lines += ["", f"Lowest common-sample WAPE: **{best_wape['model']}**. "
              f"Lowest common-sample RMSE: **{best_rmse['model']}**.", "",
              "## Information-criterion comparison groups", "",
              "DIC and WAIC came from separate full-data descriptive fits. Compare these criteria "
              "only within groups having identical fitted response observations.", ""]
    for group in dict.fromkeys(r["criteria_group"] for r in rows):
        members = [r for r in rows if r["criteria_group"] == group]
        names = ", ".join(r["model"] for r in members)
        lines.append(f"- {group}: {names}; {members[0]['fit_n']} fitted observations.")
    lines += ["", "## Sample exclusions", "", "| Model | Eligible January rows | Shared rows | Excluded |",
              "|---|---:|---:|---:|"]
    for row in audit:
        lines.append(f"| {row['model']} | {row['eligible_test_n']} | {row['common_test_n']} | {row['excluded_from_common']} |")
    lines += ["", "## Manuscript implications", "",
              "- Replace all predictive metrics and rankings using the table above.",
              "- Describe a four-week January 2023 evaluation, not a full-year or two-year holdout.",
              "- Do not compare DIC/WAIC across different fitted-response groups.",
              "- Keep full-data rainfall-effect maps and fitted-value figures distinct from held-out prediction.",
              "- Models using contemporaneous weather condition on observed weather; this is not an operational weather forecast evaluation.",
              "- Broader seasonal and epidemic generalization requires a separate rolling-origin evaluation.", ""]
    path = OUTPUT / "common_evaluation_report.md"
    path.write_text("\n".join(lines))
    print(path)

if __name__ == "__main__":
    main()
