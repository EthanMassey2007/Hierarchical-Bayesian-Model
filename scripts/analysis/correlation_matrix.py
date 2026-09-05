#!/usr/bin/env python3
"""
Create a correlation matrix for the final non-spatial model covariates.

This script does not fit a Bayesian model. It reads the combined dengue dataset,
recreates the final M5-style engineered covariates, and writes one correlation
table/figure for the variables used in the final non-spatial baseline.

Default output:
  outputs/correlation_matrix_final_model_pearson.csv
  outputs/correlation_matrix_final_model_pearson_heatmap.png
  outputs/correlation_matrix_final_model_metadata.csv

Run from the project root:
  python correlation_matrix.py
"""

from __future__ import annotations

import argparse
import math
import os
import re
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_START_YEAR = 2017
DEFAULT_END_YEAR = 2023
DEFAULT_WEATHER_LAG_WEEKS = 12
DEFAULT_CASE_LAG_WEEKS = 4

CORRELATION_VARIABLES = [
    "rainfall_lag",
    "humidity_lag",
    "temperature_lag",
    "idhm",
    "log_cases_lag",
]


def find_project_dir(project_dir_arg: str | None = None) -> Path:
    if project_dir_arg:
        return Path(project_dir_arg).expanduser().resolve()

    override = os.environ.get("HBM_PROJECT_DIR")
    if override:
        return Path(override).expanduser().resolve()

    here = Path(__file__).resolve().parent
    candidates = [Path.cwd(), here, *here.parents]
    for candidate in candidates:
        if (candidate / "data" / "complete_combined_datasets.csv").exists():
            return candidate.resolve()

    raise FileNotFoundError(
        "Could not find project directory. Run from the project root or set "
        "HBM_PROJECT_DIR=/path/to/Hierarchical-Bayesian-Model."
    )


def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [
        re.sub(r"[^0-9a-zA-Z]+", "_", str(col).strip().lower()).strip("_")
        for col in df.columns
    ]
    return df


def normalize_name(value: object) -> str:
    text = "" if pd.isna(value) else str(value)
    text = unicodedata.normalize("NFKD", text)
    text = "".join(ch for ch in text if not unicodedata.combining(ch))
    text = re.sub(r"/[a-z]{2}$", "", text.lower().strip())
    text = re.sub(r"[^0-9a-zA-Z]+", " ", text).strip()
    return re.sub(r"\s+", " ", text)


def iso_week_to_date(year: pd.Series, week: pd.Series) -> pd.Series:
    iso = (
        year.astype("Int64").astype(str)
        + "-W"
        + week.astype("Int64").astype(str).str.zfill(2)
        + "-1"
    )
    return pd.to_datetime(iso, format="%G-W%V-%u", errors="coerce")


def load_combined_data(data_dir: Path, start_year: int, end_year: int) -> pd.DataFrame:
    combined_file = data_dir / "complete_combined_datasets.csv"
    df = clean_columns(pd.read_csv(combined_file))
    required = {"municipio", "year", "week", "cases", "rainfall", "humidity", "temperature", "idhm"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Missing required columns in {combined_file}: {missing}")

    df["municipio"] = df["municipio"].map(normalize_name)
    for col in ["year", "week", "cases", "rainfall", "humidity", "temperature", "idhm"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["municipio", "year", "week", "cases"]).copy()
    df["year"] = df["year"].astype(int)
    df["week"] = df["week"].astype(int)
    df = df[(df["year"] >= start_year) & (df["year"] <= end_year)].copy()
    df["date"] = iso_week_to_date(df["year"], df["week"])
    df = df.dropna(subset=["date"]).copy()
    return df


def add_basic_lags(df: pd.DataFrame, weather_lag_weeks: int, case_lag_weeks: int) -> pd.DataFrame:
    df = df.sort_values(["municipio", "date"]).copy()

    weather_cols = ["rainfall", "humidity", "temperature"]
    weather_lookup = df[["municipio", "date", *weather_cols]].copy()
    weather_lookup["date"] += pd.to_timedelta(weather_lag_weeks * 7, unit="D")
    weather_lookup = weather_lookup.rename(columns={col: f"{col}_lag" for col in weather_cols})
    df = df.merge(weather_lookup, on=["municipio", "date"], how="left")

    case_lookup = df[["municipio", "date", "cases"]].copy()
    case_lookup["date"] += pd.to_timedelta(case_lag_weeks * 7, unit="D")
    case_lookup = case_lookup.rename(columns={"cases": "cases_lag"})
    df = df.merge(case_lookup, on=["municipio", "date"], how="left")
    df["log_cases_lag"] = np.log1p(df["cases_lag"])
    return df


def build_correlation_dataframe(args: argparse.Namespace) -> pd.DataFrame:
    project_dir = find_project_dir(args.project_dir)
    data_dir = project_dir / "data"

    df = load_combined_data(data_dir, args.start_year, args.end_year)
    df = add_basic_lags(df, args.weather_lag_weeks, args.case_lag_weeks)
    return df


def save_heatmap(corr: pd.DataFrame, output_path: Path, title: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    values = corr.to_numpy(dtype=float)
    size = max(12, 0.7 * len(corr.columns))
    fig, ax = plt.subplots(figsize=(size, size))
    image = ax.imshow(values, vmin=-1, vmax=1, cmap="coolwarm")
    ax.set_xticks(np.arange(len(corr.columns)))
    ax.set_yticks(np.arange(len(corr.index)))
    ax.set_xticklabels(corr.columns, rotation=45, ha="right", fontsize=17)
    ax.set_yticklabels(corr.index, fontsize=17)
    ax.set_title(title, fontsize = 20)

    for row_idx in range(values.shape[0]):
        for col_idx in range(values.shape[1]):
            if not np.isnan(values[row_idx, col_idx]):
                color = "white" if abs(values[row_idx, col_idx]) >= 0.55 else "black"
                ax.text(col_idx, row_idx, f"{values[row_idx, col_idx]:.2f}", ha="center", va="center", fontsize=20, color=color)

    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Correlation")
    plt.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create a final-model dengue covariate correlation matrix.")
    parser.add_argument("--project-dir", default=None)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--start-year", type=int, default=DEFAULT_START_YEAR)
    parser.add_argument("--end-year", type=int, default=DEFAULT_END_YEAR)
    parser.add_argument("--weather-lag-weeks", type=int, default=DEFAULT_WEATHER_LAG_WEEKS)
    parser.add_argument("--case-lag-weeks", type=int, default=DEFAULT_CASE_LAG_WEEKS)
    parser.add_argument("--method", choices=["pearson", "spearman"], default="pearson")
    parser.add_argument("--no-heatmap", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project_dir = find_project_dir(args.project_dir)
    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else project_dir / "outputs"
    output_dir.mkdir(parents=True, exist_ok=True)

    df = build_correlation_dataframe(args)
    variables = [col for col in CORRELATION_VARIABLES if col in df.columns]
    matrix_df = df[variables].copy()
    complete_df = matrix_df.dropna()
    if complete_df.empty:
        raise ValueError("No complete rows are available for the full correlation matrix.")

    corr = complete_df.corr(method=args.method)
    csv_path = output_dir / f"correlation_matrix_final_model_{args.method}.csv"
    png_path = output_dir / f"correlation_matrix_final_model_{args.method}_heatmap.png"
    metadata_path = output_dir / "correlation_matrix_final_model_metadata.csv"

    corr.to_csv(csv_path)
    metadata = pd.DataFrame(
        {
            "variable": variables,
            "non_missing_rows": [int(matrix_df[col].notna().sum()) for col in variables],
            "missing_rows": [int(matrix_df[col].isna().sum()) for col in variables],
            "sd": [float(matrix_df[col].std(skipna=True)) for col in variables],
        }
    )
    metadata.to_csv(metadata_path, index=False)

    if not args.no_heatmap:
        save_heatmap(corr, png_path, f"Final Model Covariates: {args.method.title()} Correlation")

    print(f"Project directory: {project_dir}")
    print(f"Rows after feature construction: {len(df)}")
    print(f"Complete rows used in matrix: {len(complete_df)}")
    print(f"Variables included: {', '.join(variables)}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {metadata_path}")
    if not args.no_heatmap:
        print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
