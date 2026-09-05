#!/usr/bin/env python3
"""
Create observed-versus-predicted dengue incidence figures for S6 and S11.

The figure compares the two candidate final models:
  S6  = best held-out WAPE/RMSE in the current model comparison.
  S11 = best WAIC/inferential fit in the current model comparison.

Full-period fitted-prediction CSVs are case-scale, so this script merges municipality population
from data/complete_combined_datasets.csv and converts observed and predicted
weekly values to incidence per 100,000 residents.

Default outputs:
  outputs/observed_predicted_figures/observed_predicted_s6_s11_data.csv
  outputs/observed_predicted_figures/observed_predicted_s6_s11.png
  outputs/observed_predicted_figures/observed_predicted_s6_s11.pdf
  outputs/observed_predicted_figures/observed_predicted_s6_s11.tiff

Run from the project root:
  python plot_observed_predicted_s6_s11.py
"""

from __future__ import annotations

import argparse
import os
import re
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd


MODEL_FILES = {
    "S6": "spatial_inla_s6_rainfall_region_full_predictions.csv",
    "S11": "s11_rainfall_spacetime_full_predictions.csv",
}

MODEL_LABELS = {
    "S6": "S6 rainfall-region",
    "S11": "S11 rainfall space-time",
}

MODEL_COLORS = {
    "Observed": "#252525",
    "S6": "#5f8f7a",
    "S11": "#9a4f5c",
}


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


def load_population(data_file: Path) -> pd.DataFrame:
    df = clean_columns(pd.read_csv(data_file, usecols=["municipio", "population"]))
    df["municipio_norm"] = df["municipio"].map(normalize_name)
    df["population"] = pd.to_numeric(df["population"], errors="coerce")
    population = (
        df.dropna(subset=["municipio_norm", "population"])
        .sort_values(["municipio_norm", "population"])
        .drop_duplicates(subset=["municipio_norm"], keep="last")
        [["municipio_norm", "population"]]
    )
    if population.empty:
        raise ValueError(f"No municipality population values found in {data_file}.")
    return population


def load_model_predictions(project_dir: Path, model: str, population: pd.DataFrame) -> pd.DataFrame:
    pred_file = project_dir / "outputs" / MODEL_FILES[model]
    if not pred_file.exists():
        raise FileNotFoundError(
            f"Missing {pred_file}. Re-run the corresponding model script to create "
            "row-level full-period fitted predictions."
        )

    df = clean_columns(pd.read_csv(pred_file))
    required = {"municipio", "year", "week", "date", "cases", "predicted_cases"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Missing required columns in {pred_file}: {missing}")

    df["municipio_norm"] = df["municipio"].map(normalize_name)
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    for col in ["year", "week", "cases", "predicted_cases"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["municipio_norm", "date", "cases", "predicted_cases"]).copy()
    df = df.merge(population, on="municipio_norm", how="left")
    missing_population = sorted(df.loc[df["population"].isna(), "municipio"].dropna().unique())
    if missing_population:
        raise ValueError(f"Missing population values for: {missing_population}")

    df["model"] = model
    df["model_label"] = MODEL_LABELS[model]
    df["observed_incidence_per_100k"] = df["cases"] / df["population"] * 100_000
    df["predicted_incidence_per_100k"] = df["predicted_cases"] / df["population"] * 100_000
    return df


def load_predictions(project_dir: Path) -> pd.DataFrame:
    population = load_population(project_dir / "data" / "complete_combined_datasets.csv")
    frames = [load_model_predictions(project_dir, model, population) for model in MODEL_FILES]
    return pd.concat(frames, ignore_index=True)


def select_representative_municipalities(predictions: pd.DataFrame, n: int = 3) -> list[str]:
    observed = (
        predictions.drop_duplicates(subset=["municipio_norm", "date"])
        .groupby("municipio_norm", as_index=False)
        .agg(total_cases=("cases", "sum"), peak_cases=("cases", "max"))
        .sort_values(["total_cases", "peak_cases"], ascending=False)
        .reset_index(drop=True)
    )
    if observed.empty:
        return []

    candidate_indices = np.linspace(0, len(observed) - 1, num=min(n, len(observed))).round().astype(int)
    return observed.iloc[candidate_indices]["municipio_norm"].tolist()


def summarize_for_plot(predictions: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    statewide = (
        predictions.groupby(["model", "model_label", "date"], as_index=False)
        .agg(
            observed_cases=("cases", "sum"),
            predicted_cases=("predicted_cases", "sum"),
            population=("population", "sum"),
        )
        .sort_values(["model", "date"])
    )
    statewide["observed_incidence_per_100k"] = (
        statewide["observed_cases"] / statewide["population"] * 100_000
    )
    statewide["predicted_incidence_per_100k"] = (
        statewide["predicted_cases"] / statewide["population"] * 100_000
    )

    scatter = (
        predictions.groupby(["model", "model_label", "municipio_norm"], as_index=False)
        .agg(
            observed_cases=("cases", "sum"),
            predicted_cases=("predicted_cases", "sum"),
            population=("population", "first"),
        )
    )
    scatter["observed_incidence_per_100k"] = scatter["observed_cases"] / scatter["population"] * 100_000
    scatter["predicted_incidence_per_100k"] = scatter["predicted_cases"] / scatter["population"] * 100_000

    selected = select_representative_municipalities(predictions)
    representative = predictions[predictions["municipio_norm"].isin(selected)].copy()
    representative["municipio_display"] = representative["municipio_norm"].str.title()
    return statewide, scatter, representative


def configure_matplotlib(output_dir: Path) -> None:
    cache_dir = output_dir / ".plot-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 600,
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def style_axis(ax) -> None:
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.6)
    ax.grid(axis="x", visible=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#4a4a4a")
    ax.spines["bottom"].set_color("#4a4a4a")
    ax.tick_params(axis="both", colors="#333333", length=3)


def save_all_formats(fig, output_stem: Path) -> None:
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_stem.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".tiff"), bbox_inches="tight")


def plot_observed_predicted(
    predictions: pd.DataFrame,
    output_dir: Path,
    show_title: bool,
) -> None:
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    statewide, scatter, representative = summarize_for_plot(predictions)

    fig = plt.figure(figsize=(9.2, 7.4))
    fig.patch.set_facecolor("white")
    grid = fig.add_gridspec(3, 2, height_ratios=[1.2, 1.0, 1.05], hspace=0.62, wspace=0.32)
    ax_time = fig.add_subplot(grid[0, :])
    ax_s6 = fig.add_subplot(grid[1, 0])
    ax_s11 = fig.add_subplot(grid[1, 1])
    representative_municipalities = representative["municipio_norm"].drop_duplicates().tolist()
    rep_grid = grid[2, :].subgridspec(
        1,
        max(1, len(representative_municipalities)),
        wspace=0.28,
    )
    rep_axes = [
        fig.add_subplot(rep_grid[0, index])
        for index in range(max(1, len(representative_municipalities)))
    ]

    if show_title:
        fig.suptitle("Observed and Predicted Dengue Incidence", x=0.02, y=0.995, ha="left", fontsize=12)

    observed_once = statewide[statewide["model"] == "S6"].copy()
    ax_time.plot(
        observed_once["date"],
        observed_once["observed_incidence_per_100k"],
        color=MODEL_COLORS["Observed"],
        linewidth=1.5,
        label="Observed",
        zorder=3,
    )
    for model in ["S6", "S11"]:
        model_dt = statewide[statewide["model"] == model]
        ax_time.plot(
            model_dt["date"],
            model_dt["predicted_incidence_per_100k"],
            color=MODEL_COLORS[model],
            linewidth=1.3,
            label=MODEL_LABELS[model],
            zorder=2,
        )
    ax_time.set_ylabel("Weekly incidence per 100,000")
    ax_time.xaxis.set_major_locator(mdates.YearLocator())
    ax_time.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax_time.set_xlim(statewide["date"].min(), statewide["date"].max())
    ax_time.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax_time.legend(frameon=False, loc="upper left", ncol=3, handlelength=2.2)
    style_axis(ax_time)

    for ax, model in [(ax_s6, "S6"), (ax_s11, "S11")]:
        model_scatter = scatter[scatter["model"] == model]
        max_value = float(
            np.nanmax(
                [
                    model_scatter["observed_incidence_per_100k"].max(),
                    model_scatter["predicted_incidence_per_100k"].max(),
                ]
            )
        )
        limit = max_value * 1.04 if max_value > 0 else 1.0
        ax.scatter(
            model_scatter["observed_incidence_per_100k"],
            model_scatter["predicted_incidence_per_100k"],
            s=18,
            color=MODEL_COLORS[model],
            alpha=0.82,
            edgecolor="white",
            linewidth=0.4,
        )
        ax.plot([0, limit], [0, limit], color="#4b4b4b", linewidth=0.8, linestyle="--")
        ax.set_xlim(0, limit)
        ax.set_ylim(0, limit)
        ax.set_title(MODEL_LABELS[model], loc="left", pad=5)
        ax.set_xlabel("Observed cumulative incidence per 100,000")
        ax.set_ylabel("Predicted cumulative incidence per 100,000")
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        style_axis(ax)

    if not representative.empty:
        for ax_rep, (municipio, municipio_dt) in zip(
            rep_axes,
            representative.groupby("municipio_norm", sort=False),
            strict=False,
        ):
            observed_dt = municipio_dt[municipio_dt["model"] == "S6"].sort_values("date")
            if observed_dt.empty:
                continue
            display = observed_dt["municipio_display"].iloc[0]
            ax_rep.plot(
                observed_dt["date"],
                observed_dt["observed_incidence_per_100k"],
                color="#9a9a9a",
                linewidth=1.0,
                alpha=0.9,
            )
            for model in ["S6", "S11"]:
                model_dt = municipio_dt[municipio_dt["model"] == model].sort_values("date")
                ax_rep.plot(
                    model_dt["date"],
                    model_dt["predicted_incidence_per_100k"],
                    color=MODEL_COLORS[model],
                    linewidth=1.0,
                    alpha=0.9,
                )
            ax_rep.set_title(display, loc="left", pad=5, fontsize=9)
            ax_rep.xaxis.set_major_locator(mdates.YearLocator(2))
            ax_rep.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
            ax_rep.set_xlim(statewide["date"].min(), statewide["date"].max())
            ax_rep.yaxis.set_major_locator(MaxNLocator(nbins=4))
            style_axis(ax_rep)

        for index, ax_rep in enumerate(rep_axes):
            if index == 0:
                ax_rep.set_ylabel("Weekly incidence per 100,000")
            else:
                ax_rep.set_ylabel("")
    else:
        rep_axes[0].axis("off")

    fig.subplots_adjust(left=0.08, right=0.98, top=0.96, bottom=0.08)
    save_all_formats(fig, output_dir / "observed_predicted_s6_s11")
    plt.close(fig)

    statewide.to_csv(output_dir / "observed_predicted_s6_s11_statewide.csv", index=False)
    scatter.to_csv(output_dir / "observed_predicted_s6_s11_scatter.csv", index=False)
    representative.to_csv(output_dir / "observed_predicted_s6_s11_representative.csv", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create observed-versus-predicted incidence figures for S6 and S11."
    )
    parser.add_argument("--project-dir", default=None)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--show-title", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project_dir = find_project_dir(args.project_dir)
    output_dir = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else project_dir / "outputs" / "observed_predicted_figures"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    configure_matplotlib(output_dir)
    predictions = load_predictions(project_dir)
    predictions.to_csv(output_dir / "observed_predicted_s6_s11_data.csv", index=False)
    plot_observed_predicted(predictions, output_dir, args.show_title)

    print(f"Rows plotted: {len(predictions):,}")
    print(f"Models plotted: {', '.join(MODEL_FILES)}")
    print(f"Wrote figures to: {output_dir}")


if __name__ == "__main__":
    main()
