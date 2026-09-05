#!/usr/bin/env python3
"""
Create a publication-ready statewide weekly dengue incidence curve.

The figure sums reported dengue cases across all municipalities in Rio de
Janeiro state for each epidemiological week in the study period, then divides by
the summed state population and multiplies by 100,000. By default it is
formatted as a study-overview figure: weekly incidence plus a 4-week rolling
mean, without model-evaluation shading.

Default outputs:
  outputs/descriptive_figures/weekly_dengue_incidence_rj_state.csv
  outputs/descriptive_figures/weekly_dengue_incidence_rj_state.png
  outputs/descriptive_figures/weekly_dengue_incidence_rj_state.pdf
  outputs/descriptive_figures/weekly_dengue_incidence_rj_state.tiff

Run from the project root:
  python plot_weekly_dengue_cases.py

To add model-evaluation shading:
  python plot_weekly_dengue_cases.py --show-test-shading
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_START_YEAR = 2017
DEFAULT_END_YEAR = 2023
DEFAULT_TEST_START = "2023-01-01"
DEFAULT_TEST_END = "2023-12-31"


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


def iso_week_to_date(year: pd.Series, week: pd.Series) -> pd.Series:
    iso = (
        year.astype("Int64").astype(str)
        + "-W"
        + week.astype("Int64").astype(str).str.zfill(2)
        + "-1"
    )
    return pd.to_datetime(iso, format="%G-W%V-%u", errors="coerce")


def load_weekly_incidence(data_file: Path, start_year: int, end_year: int) -> pd.DataFrame:
    df = clean_columns(pd.read_csv(data_file))
    required = {"municipio", "year", "week", "cases", "population"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Missing required columns in {data_file}: {missing}")

    for col in ["year", "week", "cases", "population"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["municipio", "year", "week", "cases", "population"]).copy()
    df["year"] = df["year"].astype(int)
    df["week"] = df["week"].astype(int)
    df = df[(df["year"] >= start_year) & (df["year"] <= end_year)].copy()
    df["date"] = iso_week_to_date(df["year"], df["week"])
    df = df.dropna(subset=["date"]).copy()
    df["cases"] = np.clip(np.rint(df["cases"]), 0, None).astype(int)

    population_lookup = (
        df[["municipio", "population"]]
        .drop_duplicates(subset=["municipio"])
        .assign(population=lambda x: pd.to_numeric(x["population"], errors="coerce"))
    )
    state_population = float(population_lookup["population"].sum())
    if not np.isfinite(state_population) or state_population <= 0:
        raise ValueError("State population must be positive to compute incidence.")

    weekly = (
        df.groupby(["date", "year", "week"], as_index=False)
        .agg(
            total_cases=("cases", "sum"),
            municipalities=("municipio", "nunique"),
        )
        .sort_values("date")
        .reset_index(drop=True)
    )
    weekly["state_population"] = state_population
    weekly["incidence_per_100k"] = weekly["total_cases"] / state_population * 100_000
    weekly["rolling_4wk_incidence_per_100k"] = (
        weekly["incidence_per_100k"].rolling(window=4, min_periods=1, center=True).mean()
    )
    return weekly


def save_figure(
    weekly: pd.DataFrame,
    output_png: Path,
    output_pdf: Path,
    output_tiff: Path,
    start_year: int,
    end_year: int,
    test_start: str | None,
    test_end: str | None,
    show_title: bool,
) -> None:
    cache_dir = output_png.parent / ".plot-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 600,
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, ax = plt.subplots(figsize=(7.2, 3.8), constrained_layout=True)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    if test_start and test_end:
        ax.axvspan(
            pd.to_datetime(test_start),
            pd.to_datetime(test_end),
            color="#e8e8e8",
            alpha=0.8,
            linewidth=0,
            label="Held-out evaluation period",
            zorder=0,
        )

    ax.plot(
        weekly["date"],
        weekly["incidence_per_100k"],
        color="#9a9a9a",
        linewidth=0.75,
        alpha=0.7,
        label="Weekly incidence",
        zorder=2,
    )
    ax.plot(
        weekly["date"],
        weekly["rolling_4wk_incidence_per_100k"],
        color="#8b1e2d",
        linewidth=2.0,
        label="4-week rolling mean",
        zorder=3,
    )

    if show_title:
        ax.set_title(f"Weekly Dengue Incidence in Rio de Janeiro State, {start_year}-{end_year}", loc="left", pad=8)
    ax.set_xlabel("Year")
    ax.set_ylabel("Weekly dengue incidence per 100,000 residents")
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.1f}"))
    year_ticks = [pd.Timestamp(year=year, month=1, day=1) for year in range(start_year, end_year + 1)]
    ax.set_xticks(year_ticks)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[4, 7, 10]))
    ax.set_xlim(pd.Timestamp(year=start_year, month=1, day=1), pd.Timestamp(year=end_year, month=12, day=31))

    ax.grid(axis="y", color="#d9d9d9", linewidth=0.6)
    ax.grid(axis="x", visible=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#4a4a4a")
    ax.spines["bottom"].set_color("#4a4a4a")
    ax.tick_params(axis="both", colors="#333333", length=3)
    ax.margins(x=0.01)

    legend = ax.legend(
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0.0, 0.98),
        handlelength=2.4,
    )
    for line in legend.get_lines():
        line.set_linewidth(2.0)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, bbox_inches="tight")
    fig.savefig(output_pdf, bbox_inches="tight")
    fig.savefig(output_tiff, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a publication-ready statewide dengue incidence time-series figure."
    )
    parser.add_argument("--project-dir", default=None)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--start-year", type=int, default=DEFAULT_START_YEAR)
    parser.add_argument("--end-year", type=int, default=DEFAULT_END_YEAR)
    parser.add_argument("--test-start", default=DEFAULT_TEST_START)
    parser.add_argument("--test-end", default=DEFAULT_TEST_END)
    parser.add_argument("--show-test-shading", action="store_true")
    parser.add_argument("--show-title", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project_dir = find_project_dir(args.project_dir)
    data_file = project_dir / "data" / "complete_combined_datasets.csv"
    output_dir = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else project_dir / "outputs" / "descriptive_figures"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    weekly = load_weekly_incidence(data_file, args.start_year, args.end_year)
    csv_path = output_dir / "weekly_dengue_incidence_rj_state.csv"
    png_path = output_dir / "weekly_dengue_incidence_rj_state.png"
    pdf_path = output_dir / "weekly_dengue_incidence_rj_state.pdf"
    tiff_path = output_dir / "weekly_dengue_incidence_rj_state.tiff"

    weekly.to_csv(csv_path, index=False)
    save_figure(
        weekly,
        png_path,
        pdf_path,
        tiff_path,
        args.start_year,
        args.end_year,
        args.test_start if args.show_test_shading else None,
        args.test_end if args.show_test_shading else None,
        args.show_title,
    )

    print(f"Project directory: {project_dir}")
    print(f"Weeks plotted: {len(weekly)}")
    print(f"Municipalities per week: {int(weekly['municipalities'].min())}-{int(weekly['municipalities'].max())}")
    print(f"State population denominator: {int(weekly['state_population'].iloc[0]):,}")
    print(f"Total cases: {int(weekly['total_cases'].sum()):,}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {png_path}")
    print(f"Wrote {pdf_path}")
    print(f"Wrote {tiff_path}")


if __name__ == "__main__":
    main()
