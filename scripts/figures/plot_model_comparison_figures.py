#!/usr/bin/env python3
"""
Create publication-ready model-comparison figures for M0-M5 and S1-S11.

The script reads outputs/all_model_results_table.csv, standardizes model labels,
adds stage groupings along the x-axis, and writes PNG, PDF, TIFF, and CSV files
to a dedicated model-comparison figure folder.

Default outputs:
  outputs/model_comparison_figures/model_comparison_data.csv
  outputs/model_comparison_figures/model_comparison_primary.png
  outputs/model_comparison_figures/model_comparison_primary.pdf
  outputs/model_comparison_figures/model_comparison_primary.tiff
  outputs/model_comparison_figures/model_comparison_heldout_metrics.png
  outputs/model_comparison_figures/model_comparison_heldout_metrics.pdf
  outputs/model_comparison_figures/model_comparison_heldout_metrics.tiff
  outputs/model_comparison_figures/model_comparison_information_criteria.png
  outputs/model_comparison_figures/model_comparison_information_criteria.pdf
  outputs/model_comparison_figures/model_comparison_information_criteria.tiff

Run from the project root:
  python plot_model_comparison_figures.py
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


MODEL_ORDER = [
    "M0",
    "M1",
    "M2",
    "M3",
    "M4",
    "M5",
    "S1",
    "S2",
    "S3",
    "S4",
    "S5",
    "S6",
    "S7",
    "S8",
    "S9",
    "S10",
    "S11",
]

STAGES = [
    ("Climate/socioeconomic", ["M0", "M1", "M2", "M3"]),
    ("Temporal", ["M4", "M5"]),
    ("Spatial", ["S1", "S2", "S3"]),
    ("Mobility", ["S4", "S5"]),
    ("Climate heterogeneity", ["S6", "S7", "S8", "S9", "S10", "S11"]),
]

STAGE_COLORS = {
    "Climate/socioeconomic": "#6f7f8f",
    "Temporal": "#b07b3f",
    "Spatial": "#5f8f7a",
    "Mobility": "#7b6fa6",
    "Climate heterogeneity": "#9a4f5c",
}

STAGE_DISPLAY_LABELS = {
    "Climate/socioeconomic": "Climate and\nsocioeconomic",
    "Temporal": "Temporal",
    "Spatial": "Spatial",
    "Mobility": "Mobility",
    "Climate heterogeneity": "Climate\nheterogeneity",
}

METRIC_SPECS = [
    ("mae", "Held-out MAE", "Mean absolute error", False),
    ("rmse", "Held-out RMSE", "Root mean squared error", False),
    ("wape_percent", "Held-out WAPE (%)", "Weighted absolute percentage error (%)", False),
    ("r2", "Held-out R$^2$", "Coefficient of determination", True),
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
        if (candidate / "outputs" / "all_model_results_table.csv").exists():
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


def clean_model_label(value: object) -> str:
    label = str(value).strip()
    if label.startswith("R_M"):
        label = label.replace("R_M", "M", 1)
    return label


def stage_for_model(model: str) -> str:
    for stage, models in STAGES:
        if model in models:
            return stage
    raise ValueError(f"Model {model!r} is not assigned to a stage.")


def load_model_results(results_file: Path) -> pd.DataFrame:
    df = clean_columns(pd.read_csv(results_file))
    required = {"model", "dic", "waic", "mae", "rmse", "wape", "r2",
                "evaluation_protocol", "evaluation_signature", "test_n", "test_start", "test_end",
                "fit_signature", "fit_n", "criteria_group", "run_id"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Missing required columns in {results_file}: {missing}")

    df["model_label"] = df["model"].map(clean_model_label)
    unexpected = sorted(set(df["model_label"]).difference(MODEL_ORDER))
    if unexpected:
        raise ValueError(f"Unexpected model labels in {results_file}: {unexpected}")

    numeric_cols = ["dic", "waic", "mae", "rmse", "wape", "r2"]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if df[["model_label", *numeric_cols]].isna().any().any() or not np.isfinite(df[numeric_cols].to_numpy()).all():
        raise ValueError("Incomplete model results; refusing to silently drop models.")
    if df["model_label"].duplicated().any() or set(df["model_label"]) != set(MODEL_ORDER):
        raise ValueError("Expected exactly one row for each of the 17 models.")
    for col in ["evaluation_protocol", "evaluation_signature", "test_n", "test_start", "test_end", "run_id"]:
        if df[col].isna().any() or df[col].nunique() != 1:
            raise ValueError(f"Models do not share the same evaluation contract: {col}")
    if df["evaluation_protocol"].iloc[0] != "common_observations_2023_january_v1":
        raise ValueError("Unrecognized evaluation protocol.")

    order_lookup = {model: index for index, model in enumerate(MODEL_ORDER)}
    df["model_index"] = df["model_label"].map(order_lookup)
    df = df.sort_values("model_index").reset_index(drop=True)
    df["stage"] = df["model_label"].map(stage_for_model)
    df["wape_percent"] = df["wape"] * 100
    df["delta_waic"] = df["waic"] - df.groupby("fit_signature")["waic"].transform("min")
    df["delta_dic"] = df["dic"] - df.groupby("fit_signature")["dic"].transform("min")
    df["waic_rank_recomputed"] = df.groupby("fit_signature")["waic"].rank(method="min").astype(int)
    df["wape_rank_recomputed"] = df["wape"].rank(method="min").astype(int)
    df["rmse_rank_recomputed"] = df["rmse"].rank(method="min").astype(int)
    return df


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


def model_colors(models: list[str]) -> list[str]:
    return [STAGE_COLORS[stage_for_model(model)] for model in models]


def add_stage_guides(ax, models: list[str], y_text: float = -0.29, y_line: float = -0.18) -> None:
    transform = ax.get_xaxis_transform()
    for stage, stage_models in STAGES:
        positions = [models.index(model) for model in stage_models if model in models]
        if not positions:
            continue
        start = min(positions) - 0.42
        end = max(positions) + 0.42
        center = (start + end) / 2
        color = STAGE_COLORS[stage]
        ax.plot(
            [start, end],
            [y_line, y_line],
            transform=transform,
            color=color,
            linewidth=1.8,
            solid_capstyle="butt",
            clip_on=False,
        )
        ax.text(
            center,
            y_text,
            STAGE_DISPLAY_LABELS[stage],
            transform=transform,
            ha="center",
            va="top",
            color="#333333",
            fontsize=7.0,
            linespacing=1.05,
            clip_on=False,
        )


def style_model_axis(ax, models: list[str], ylabel: str) -> None:
    ax.set_xlim(-0.6, len(models) - 0.4)
    ax.set_xticks(range(len(models)))
    ax.set_xticklabels(models, fontsize=8)
    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.6)
    ax.grid(axis="x", visible=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#4a4a4a")
    ax.spines["bottom"].set_color("#4a4a4a")
    ax.tick_params(axis="both", colors="#333333", length=3)
    add_stage_guides(ax, models)


def save_all_formats(fig, output_stem: Path) -> None:
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_stem.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".tiff"), bbox_inches="tight")


def plot_line_metric(ax, df: pd.DataFrame, metric: str, ylabel: str, higher_is_better: bool = False) -> None:
    models = df["model_label"].tolist()
    x = np.arange(len(df))
    colors = model_colors(models)
    y = df[metric].to_numpy(dtype=float)
    best_idx = int(np.argmax(y) if higher_is_better else np.argmin(y))

    ax.plot(x, y, color="#4b4b4b", linewidth=0.9, zorder=1)
    ax.scatter(x, y, s=34, c=colors, edgecolor="white", linewidth=0.8, zorder=2)
    ax.scatter(
        [x[best_idx]],
        [y[best_idx]],
        s=72,
        facecolor="none",
        edgecolor="#111111",
        linewidth=1.2,
        zorder=3,
    )
    style_model_axis(ax, models, ylabel)


def plot_primary(df: pd.DataFrame, output_dir: Path, show_title: bool) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(8.2, 6.2), sharex=False)
    fig.patch.set_facecolor("white")
    if show_title:
        fig.suptitle("Model Comparison", x=0.02, y=0.99, ha="left", fontsize=12, fontweight="bold")

    plot_line_metric(axes[0], df, "rmse", "Held-out RMSE", higher_is_better=False)
    axes[0].set_xlabel("")
    axes[0].set_xticklabels([])
    for text in list(axes[0].texts):
        text.remove()
    add_stage_guides(axes[0], df["model_label"].tolist(), y_text=-0.24, y_line=-0.13)

    plot_line_metric(axes[1], df, "wape_percent", "Held-out WAPE (%)", higher_is_better=False)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.95, bottom=0.24, hspace=0.65)
    save_all_formats(fig, output_dir / "model_comparison_primary")
    plt.close(fig)


def plot_heldout_metrics(df: pd.DataFrame, output_dir: Path, show_title: bool) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.8))
    fig.patch.set_facecolor("white")
    if show_title:
        fig.suptitle("Held-out Predictive Performance", x=0.02, y=0.99, ha="left", fontsize=12, fontweight="bold")

    for ax, (metric, _, ylabel, higher_is_better) in zip(axes.flat, METRIC_SPECS, strict=True):
        plot_line_metric(ax, df, metric, ylabel, higher_is_better=higher_is_better)

    fig.text(0.5, 0.02, f"Common sample: n={int(df['test_n'].iloc[0])}; "
             f"{df['test_start'].iloc[0]} to {df['test_end'].iloc[0]}", ha="center", fontsize=9)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.20, hspace=0.85, wspace=0.28)
    save_all_formats(fig, output_dir / "model_comparison_heldout_metrics")
    plt.close(fig)


def plot_information_criteria(df: pd.DataFrame, output_dir: Path, show_title: bool) -> None:
    import matplotlib.pyplot as plt

    groups = [group for _, group in df.groupby("fit_signature", sort=False) if len(group) > 1]
    if not groups:
        raise ValueError("No models share fitted observations; DIC/WAIC cannot be compared.")
    fig, axes = plt.subplots(len(groups), 1, figsize=(8.2, 4.2 * len(groups)), squeeze=False)
    for ax, group in zip(axes.flat, groups, strict=True):
        models = group["model_label"].tolist()
        x = np.arange(len(group))
        ax.plot(x, group["delta_waic"], marker="o", label="Delta WAIC")
        ax.plot(x, group["delta_dic"], marker="o", linestyle="--", label="Delta DIC")
        ax.set_xticks(x, models)
        ax.set_ylabel("Difference within fitted sample")
        ax.set_title(f"{int(group['fit_n'].iloc[0]):,} identical fitted observations")
        ax.legend(frameon=False)
    fig.tight_layout()
    save_all_formats(fig, output_dir / "model_comparison_information_criteria")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create grouped, publication-ready model-comparison figures."
    )
    parser.add_argument("--project-dir", default=None)
    parser.add_argument("--results-file", default=None)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--show-title", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project_dir = find_project_dir(args.project_dir)
    results_file = (
        Path(args.results_file).expanduser().resolve()
        if args.results_file
        else project_dir / "outputs" / "all_model_results_table.csv"
    )
    output_dir = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else project_dir / "outputs" / "model_comparison_figures"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    configure_matplotlib(output_dir)
    status_file = project_dir / "outputs" / "evaluation_status.csv"
    if not status_file.exists():
        raise ValueError("No verified evaluation run. Run make all-results first.")
    status = pd.read_csv(status_file)
    if len(status) != 1 or status["status"].iloc[0] != "complete":
        raise ValueError("Evaluation is incomplete or failed; old figures must not be regenerated.")
    results = load_model_results(results_file)
    if results["run_id"].iloc[0] != status["run_id"].iloc[0]:
        raise ValueError("Results do not match the completed evaluation run.")
    results.to_csv(output_dir / "model_comparison_data.csv", index=False)

    plot_primary(results, output_dir, args.show_title)
    plot_heldout_metrics(results, output_dir, args.show_title)
    plot_information_criteria(results, output_dir, args.show_title)

    best_wape = results.loc[results["wape"].idxmin(), "model_label"]
    best_rmse = results.loc[results["rmse"].idxmin(), "model_label"]
    print(f"Models plotted: {', '.join(results['model_label'])}")
    print("DIC/WAIC comparisons are restricted to identical fitted-response samples.")
    print(f"Best held-out WAPE: {best_wape}")
    print(f"Best held-out RMSE: {best_rmse}")
    print(f"Wrote figures to: {output_dir}")


if __name__ == "__main__":
    main()
