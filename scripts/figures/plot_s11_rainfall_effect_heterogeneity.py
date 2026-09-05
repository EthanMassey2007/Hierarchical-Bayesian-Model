#!/usr/bin/env python3
"""
Create the main rainfall-effect heterogeneity figure from the final S11 model.

The figure contains:
  (a) time-averaged municipality-level rainfall relative risk (RR)
  (b) temporal rainfall RR with an approximate 95% credible interval
  (c) municipality-by-time rainfall RR heatmap ordered by IBGE intermediate region
  (d) four representative periods showing how spatial rainfall sensitivity changes

Default outputs:
  outputs/rainfall_effect_heterogeneity_figures/s11_rainfall_effect_heterogeneity.png
  outputs/rainfall_effect_heterogeneity_figures/s11_rainfall_effect_heterogeneity.pdf
  outputs/rainfall_effect_heterogeneity_figures/s11_rainfall_effect_heterogeneity.tiff

Run from the project root:
  python plot_s11_rainfall_effect_heterogeneity.py
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_OUTPUT_NAME = "s11_rainfall_effect_heterogeneity"
REGION_FIELD = "regiao_intermediaria_nome"


def find_project_dir(project_dir_arg: str | None = None) -> Path:
    if project_dir_arg:
        return Path(project_dir_arg).expanduser().resolve()

    override = os.environ.get("HBM_PROJECT_DIR")
    if override:
        return Path(override).expanduser().resolve()

    here = Path(__file__).resolve().parent
    candidates = [Path.cwd(), here, *here.parents]
    for candidate in candidates:
        if (candidate / "outputs" / "s11_rainfall_spacetime_combined_effects.csv").exists():
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


def polygon_rings(geometry: dict) -> list[list[list[float]]]:
    geom_type = geometry.get("type")
    coordinates = geometry.get("coordinates", [])
    if geom_type == "Polygon":
        return [coordinates[0]]
    if geom_type == "MultiPolygon":
        return [polygon[0] for polygon in coordinates]
    return []


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
            "font.size": 8.5,
            "axes.titlesize": 9.5,
            "axes.labelsize": 8.5,
            "axes.linewidth": 0.75,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def load_region_lookup(data_dir: Path) -> pd.DataFrame:
    path = data_dir / "rj_ibge_intermediate_regions.csv"
    regions = clean_columns(pd.read_csv(path))
    required = {"ibge_code", "municipio", REGION_FIELD}
    missing = sorted(required.difference(regions.columns))
    if missing:
        raise ValueError(f"Missing required columns in {path}: {missing}")

    regions["ibge_code"] = pd.to_numeric(regions["ibge_code"], errors="coerce").astype("Int64")
    regions["municipio_norm"] = regions["municipio"].map(normalize_name)
    regions = regions.dropna(subset=["ibge_code", REGION_FIELD]).copy()
    regions["ibge_code"] = regions["ibge_code"].astype(int)
    regions[REGION_FIELD] = regions[REGION_FIELD].astype(str)
    return regions[["ibge_code", "municipio_norm", REGION_FIELD]].drop_duplicates("ibge_code")


def load_s11_effects(project_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    output_dir = project_dir / "outputs"
    combined_path = output_dir / "s11_rainfall_spacetime_combined_effects.csv"
    time_path = output_dir / "s11_rainfall_spacetime_time_effects.csv"
    muni_path = output_dir / "s11_rainfall_spacetime_municipality_effects.csv"
    for path in [combined_path, time_path, muni_path]:
        if not path.exists():
            raise FileNotFoundError(
                f"Missing {path}. Re-run models/r_inla/spatial/spatial_inla_model_s11_rainfall_spacetime.R."
            )

    combined = clean_columns(pd.read_csv(combined_path))
    time_effects = clean_columns(pd.read_csv(time_path))
    municipality_effects = clean_columns(pd.read_csv(muni_path))

    for name, df in {
        "combined": combined,
        "time_effects": time_effects,
        "municipality_effects": municipality_effects,
    }.items():
        required = {"ibge_code", "rainfall_relative_risk"} if name != "time_effects" else {
            "date",
            "rainfall_relative_risk",
            "rainfall_relative_risk_q025",
            "rainfall_relative_risk_q975",
        }
        missing = sorted(required.difference(df.columns))
        if missing:
            raise ValueError(f"Missing required columns in {name}: {missing}")

    combined["date"] = pd.to_datetime(combined["date"], errors="coerce")
    time_effects["date"] = pd.to_datetime(time_effects["date"], errors="coerce")
    for df in [combined, time_effects, municipality_effects]:
        for col in df.columns:
            if col.startswith("rainfall_") or col in {"ibge_code", "year", "week"}:
                df[col] = pd.to_numeric(df[col], errors="coerce")

    combined = combined.dropna(subset=["ibge_code", "date", "rainfall_relative_risk"]).copy()
    combined["ibge_code"] = combined["ibge_code"].astype(int)
    municipality_effects = municipality_effects.dropna(subset=["ibge_code", "rainfall_relative_risk"]).copy()
    municipality_effects["ibge_code"] = municipality_effects["ibge_code"].astype(int)
    time_effects = time_effects.dropna(
        subset=[
            "date",
            "rainfall_relative_risk",
            "rainfall_relative_risk_q025",
            "rainfall_relative_risk_q975",
        ]
    ).copy()

    return combined, time_effects, municipality_effects


def prepare_effect_tables(
    combined: pd.DataFrame,
    time_effects: pd.DataFrame,
    municipality_effects: pd.DataFrame,
    regions: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    combined = combined.merge(regions, on="ibge_code", how="left")
    municipality_effects = municipality_effects.merge(regions, on="ibge_code", how="left")
    missing_regions = sorted(combined.loc[combined[REGION_FIELD].isna(), "ibge_code"].unique())
    if missing_regions:
        raise ValueError(f"Missing IBGE intermediate-region lookup for codes: {missing_regions}")

    spatial_average = (
        combined.groupby(["ibge_code", "municipio", REGION_FIELD], as_index=False)
        .agg(
            rainfall_relative_risk=("rainfall_relative_risk", "mean"),
            rainfall_relative_risk_q025=("rainfall_relative_risk_q025", "mean"),
            rainfall_relative_risk_q975=("rainfall_relative_risk_q975", "mean"),
        )
        .sort_values("rainfall_relative_risk", ascending=False)
    )

    order = (
        spatial_average.assign(municipio_norm=lambda x: x["municipio"].map(normalize_name))
        .sort_values([REGION_FIELD, "rainfall_relative_risk", "municipio_norm"], ascending=[True, False, True])
        [["ibge_code", "municipio", REGION_FIELD, "rainfall_relative_risk"]]
        .reset_index(drop=True)
    )
    order["municipality_order"] = np.arange(len(order))

    heatmap = combined.merge(order[["ibge_code", "municipality_order"]], on="ibge_code", how="inner")
    heatmap = heatmap.sort_values(["municipality_order", "date"])

    time_effects = time_effects.sort_values("date").copy()
    return spatial_average, time_effects, heatmap, order


def build_map_patches(geojson: dict, values_by_ibge: dict[int, float]):
    from matplotlib.patches import Polygon

    patches = []
    values = []
    for feature in geojson["features"]:
        props = feature.get("properties", {})
        raw_code = props.get("GEOCODIGO")
        if raw_code is None or pd.isna(raw_code):
            continue
        ibge_code = int(raw_code)
        value = values_by_ibge.get(ibge_code)
        if value is None or not math.isfinite(value):
            continue
        for ring in polygon_rings(feature["geometry"]):
            patches.append(Polygon(ring, closed=True))
            values.append(value)

    if not patches:
        raise ValueError("No GeoJSON polygons matched S11 rainfall-effect values.")
    return patches, np.asarray(values, dtype=float)


def add_map(
    ax,
    geojson: dict,
    values_by_ibge: dict[int, float],
    norm,
    cmap: str,
    title: str,
    linewidth: float = 0.12,
):
    from matplotlib.collections import PatchCollection

    patches, values = build_map_patches(geojson, values_by_ibge)
    collection = PatchCollection(
        patches,
        cmap=cmap,
        norm=norm,
        edgecolor="#ffffff",
        linewidth=linewidth,
    )
    collection.set_array(values)
    ax.add_collection(collection)
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, loc="left", pad=3)
    return collection


def style_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#4a4a4a")
    ax.spines["bottom"].set_color("#4a4a4a")
    ax.tick_params(axis="both", colors="#333333", length=3)
    ax.grid(axis="y", color="#dedede", linewidth=0.55)
    ax.grid(axis="x", visible=False)


def save_all_formats(fig, output_stem: Path) -> None:
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_stem.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".tiff"), bbox_inches="tight")


def representative_dates(time_effects: pd.DataFrame, n: int = 4) -> list[pd.Timestamp]:
    dates = time_effects["date"].drop_duplicates().sort_values().reset_index(drop=True)
    if dates.empty:
        return []
    positions = np.linspace(0.12, 0.88, n)
    indices = np.clip(np.rint(positions * (len(dates) - 1)).astype(int), 0, len(dates) - 1)
    return list(dates.iloc[indices].drop_duplicates())


def plot_figure(
    project_dir: Path,
    output_dir: Path,
    spatial_average: pd.DataFrame,
    time_effects: pd.DataFrame,
    heatmap: pd.DataFrame,
    order: pd.DataFrame,
    show_title: bool,
) -> None:
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib.ticker import MaxNLocator

    with (project_dir / "data" / "RJ.json").open(encoding="utf-8") as f:
        geojson = json.load(f)

    rr_values = pd.concat(
        [
            spatial_average["rainfall_relative_risk"],
            time_effects["rainfall_relative_risk"],
            heatmap["rainfall_relative_risk"],
        ],
        ignore_index=True,
    ).dropna()
    lower = float(max(0.2, np.nanquantile(rr_values, 0.01)))
    upper = float(np.nanquantile(rr_values, 0.99))
    span = max(upper - 1.0, 1.0 - lower)
    norm = TwoSlopeNorm(vmin=max(0.05, 1.0 - span), vcenter=1.0, vmax=1.0 + span)
    cmap = "RdBu_r"

    fig = plt.figure(figsize=(10.6, 9.0))
    fig.patch.set_facecolor("white")
    grid = fig.add_gridspec(
        3,
        2,
        height_ratios=[1.0, 1.18, 0.72],
        width_ratios=[0.92, 1.08],
        hspace=0.52,
        wspace=0.23,
    )
    ax_map = fig.add_subplot(grid[0, 0])
    ax_time = fig.add_subplot(grid[0, 1])
    ax_heat = fig.add_subplot(grid[1, :])
    period_grid = grid[2, :].subgridspec(1, 4, wspace=0.08)
    period_axes = [fig.add_subplot(period_grid[0, i]) for i in range(4)]

    if show_title:
        fig.suptitle("S11 Rainfall-Effect Heterogeneity", x=0.03, y=0.995, ha="left", fontsize=12)

    map_values = dict(
        zip(
            spatial_average["ibge_code"].astype(int),
            spatial_average["rainfall_relative_risk"].astype(float),
            strict=True,
        )
    )
    map_collection = add_map(ax_map, geojson, map_values, norm, cmap, "(a) Time-averaged rainfall RR")
    cbar = fig.colorbar(map_collection, ax=ax_map, fraction=0.035, pad=0.015)
    cbar.outline.set_visible(False)
    cbar.set_label("RR per 1-SD rainfall increase", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    ax_time.fill_between(
        time_effects["date"],
        time_effects["rainfall_relative_risk_q025"],
        time_effects["rainfall_relative_risk_q975"],
        color="#b8c7d9",
        alpha=0.7,
        linewidth=0,
    )
    ax_time.plot(time_effects["date"], time_effects["rainfall_relative_risk"], color="#2b5c7c", linewidth=1.35)
    ax_time.axhline(1, color="#4a4a4a", linewidth=0.85, linestyle="--")
    ax_time.set_title("(b) Temporal rainfall RR", loc="left", pad=5)
    ax_time.set_ylabel("RR per 1-SD rainfall increase")
    ax_time.set_xlabel("Year")
    ax_time.set_xlim(time_effects["date"].min(), time_effects["date"].max())
    ax_time.xaxis.set_major_locator(mdates.YearLocator())
    ax_time.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax_time.yaxis.set_major_locator(MaxNLocator(nbins=5))
    style_axis(ax_time)

    heat_matrix = (
        heatmap.pivot_table(index="municipality_order", columns="date", values="rainfall_relative_risk", aggfunc="mean")
        .sort_index()
    )
    heat_dates = pd.to_datetime(heat_matrix.columns)
    extent = [
        mdates.date2num(heat_dates.min()),
        mdates.date2num(heat_dates.max()),
        heat_matrix.index.max() + 0.5,
        heat_matrix.index.min() - 0.5,
    ]
    image = ax_heat.imshow(
        heat_matrix.to_numpy(),
        aspect="auto",
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
        extent=extent,
    )
    ax_heat.set_title(
        "(c) Municipality x time rainfall RR, ordered by IBGE intermediate region",
        loc="left",
        pad=5,
    )
    ax_heat.set_xlabel("Year")
    ax_heat.set_ylabel("")
    ax_heat.xaxis_date()
    ax_heat.xaxis.set_major_locator(mdates.YearLocator())
    ax_heat.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax_heat.set_yticks([])
    ax_heat.spines["top"].set_visible(False)
    ax_heat.spines["right"].set_visible(False)
    ax_heat.spines["left"].set_color("#4a4a4a")
    ax_heat.spines["bottom"].set_color("#4a4a4a")
    heat_cbar = fig.colorbar(image, ax=ax_heat, fraction=0.018, pad=0.01)
    heat_cbar.outline.set_visible(False)
    heat_cbar.set_label("RR", fontsize=8)
    heat_cbar.ax.tick_params(labelsize=7)

    region_ranges = (
        order.groupby(REGION_FIELD, sort=False)["municipality_order"]
        .agg(["min", "max"])
        .reset_index()
    )
    for _, row in region_ranges.iterrows():
        boundary = row["max"] + 0.5
        ax_heat.axhline(boundary, color="#ffffff", linewidth=0.55, alpha=0.9)
        midpoint = (row["min"] + row["max"]) / 2
        ax_heat.text(
            mdates.date2num(heat_dates.min()) - (mdates.date2num(heat_dates.max()) - mdates.date2num(heat_dates.min())) * 0.012,
            midpoint,
            str(row[REGION_FIELD]),
            ha="right",
            va="center",
            fontsize=6.5,
            color="#333333",
            clip_on=False,
        )

    dates = representative_dates(time_effects, n=4)
    for ax_period, date in zip(period_axes, dates, strict=False):
        period = heatmap[heatmap["date"] == date]
        values = dict(zip(period["ibge_code"].astype(int), period["rainfall_relative_risk"].astype(float), strict=True))
        add_map(
            ax_period,
            geojson,
            values,
            norm,
            cmap,
            f"{date.strftime('%b %Y')}",
            linewidth=0.08,
        )
    for ax_period in period_axes[len(dates):]:
        ax_period.axis("off")
    period_axes[0].text(
        0,
        1.24,
        "(d) Spatial rainfall RR at representative periods",
        transform=period_axes[0].transAxes,
        ha="left",
        va="bottom",
        fontsize=9.5,
    )

    fig.subplots_adjust(left=0.13, right=0.985, top=0.965, bottom=0.055)
    save_all_formats(fig, output_dir / DEFAULT_OUTPUT_NAME)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create the S11 rainfall-effect heterogeneity figure."
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
        else project_dir / "outputs" / "rainfall_effect_heterogeneity_figures"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    configure_matplotlib(output_dir)

    combined, time_effects, municipality_effects = load_s11_effects(project_dir)
    regions = load_region_lookup(project_dir / "data")
    spatial_average, time_effects, heatmap, order = prepare_effect_tables(
        combined,
        time_effects,
        municipality_effects,
        regions,
    )

    spatial_average.to_csv(output_dir / "s11_rainfall_rr_spatial_average.csv", index=False)
    time_effects.to_csv(output_dir / "s11_rainfall_rr_temporal_curve.csv", index=False)
    heatmap.to_csv(output_dir / "s11_rainfall_rr_municipality_time_heatmap.csv", index=False)
    order.to_csv(output_dir / "s11_rainfall_rr_municipality_order.csv", index=False)

    plot_figure(project_dir, output_dir, spatial_average, time_effects, heatmap, order, args.show_title)

    print(f"Municipalities plotted: {spatial_average['ibge_code'].nunique():,}")
    print(f"Weeks plotted: {time_effects['date'].nunique():,}")
    print(f"Municipality-week effects plotted: {len(heatmap):,}")
    print(f"Wrote figures to: {output_dir}")


if __name__ == "__main__":
    main()
