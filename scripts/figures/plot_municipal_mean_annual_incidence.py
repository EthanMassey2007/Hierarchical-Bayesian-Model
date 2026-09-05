#!/usr/bin/env python3
"""
Create a publication-ready map of municipal mean annual dengue incidence.

For each Rio de Janeiro municipality, the script sums reported cases by year,
averages those annual totals across the study period, and divides by the
municipality population to compute mean annual incidence per 100,000 residents.

Default outputs:
  outputs/descriptive_figures/municipal_mean_annual_dengue_incidence.csv
  outputs/descriptive_figures/municipal_mean_annual_dengue_incidence.png
  outputs/descriptive_figures/municipal_mean_annual_dengue_incidence.pdf
  outputs/descriptive_figures/municipal_mean_annual_dengue_incidence.tiff

Run from the project root:
  python plot_municipal_mean_annual_incidence.py
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


DEFAULT_START_YEAR = 2017
DEFAULT_END_YEAR = 2023


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


def load_municipio_lookup(data_dir: Path) -> pd.DataFrame:
    region_lookup_file = data_dir / "rj_ibge_intermediate_regions.csv"
    if region_lookup_file.exists():
        regions = clean_columns(pd.read_csv(region_lookup_file))
        required = {"municipio", "ibge_code"}
        missing = sorted(required.difference(regions.columns))
        if missing:
            raise ValueError(f"Missing required columns in {region_lookup_file}: {missing}")
        regions["municipio_norm"] = regions["municipio"].map(normalize_name)
        regions["ibge_code"] = pd.to_numeric(regions["ibge_code"], errors="coerce").astype("Int64")
        regions = regions[regions["ibge_code"].notna()].drop_duplicates(subset=["municipio_norm"])
        return regions[["municipio_norm", "ibge_code"]]

    municipios_file = data_dir / "municipios.csv"
    municipios = clean_columns(pd.read_csv(municipios_file, usecols=["city", "ibgeID", "state"]))
    municipios["state_code"] = municipios["state"].astype(str).str.upper()
    municipios["municipio_norm"] = municipios["city"].map(normalize_name)
    municipios["ibge_code"] = pd.to_numeric(municipios["ibgeid"], errors="coerce").astype("Int64")
    municipios = municipios[
        (municipios["state_code"] == "RJ") & municipios["ibge_code"].notna()
    ].copy()
    municipios = municipios.drop_duplicates(subset=["municipio_norm"])
    return municipios[["municipio_norm", "ibge_code"]]


def load_mean_annual_incidence(data_dir: Path, start_year: int, end_year: int) -> pd.DataFrame:
    data_file = data_dir / "complete_combined_datasets.csv"
    df = clean_columns(pd.read_csv(data_file))
    required = {"municipio", "year", "cases", "population"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Missing required columns in {data_file}: {missing}")

    for col in ["year", "cases", "population"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["municipio", "year", "cases", "population"]).copy()
    df["year"] = df["year"].astype(int)
    df = df[(df["year"] >= start_year) & (df["year"] <= end_year)].copy()
    df["municipio_norm"] = df["municipio"].map(normalize_name)
    df["cases"] = np.clip(np.rint(df["cases"]), 0, None).astype(int)

    lookup = load_municipio_lookup(data_dir)
    df = df.merge(lookup, on="municipio_norm", how="left")
    missing_codes = sorted(df.loc[df["ibge_code"].isna(), "municipio"].dropna().unique())
    if missing_codes:
        raise ValueError(f"Municipalities missing IBGE codes: {missing_codes}")

    annual = (
        df.groupby(["ibge_code", "municipio_norm", "year"], as_index=False)
        .agg(
            annual_cases=("cases", "sum"),
            population=("population", "first"),
            weekly_rows=("cases", "size"),
        )
    )
    years_expected = end_year - start_year + 1
    summary = (
        annual.groupby(["ibge_code", "municipio_norm"], as_index=False)
        .agg(
            mean_annual_cases=("annual_cases", "mean"),
            total_cases=("annual_cases", "sum"),
            population=("population", "first"),
            years_observed=("year", "nunique"),
            mean_weekly_rows_per_year=("weekly_rows", "mean"),
        )
    )
    summary["years_expected"] = years_expected
    summary["mean_annual_incidence_per_100k"] = (
        summary["mean_annual_cases"] / summary["population"] * 100_000
    )
    summary["cumulative_incidence_per_100k"] = (
        summary["total_cases"] / summary["population"] * 100_000
    )
    summary["ibge_code"] = summary["ibge_code"].astype(int)
    return summary.sort_values("mean_annual_incidence_per_100k", ascending=False)


def polygon_rings(geometry: dict) -> list[list[list[float]]]:
    geom_type = geometry.get("type")
    coordinates = geometry.get("coordinates", [])
    if geom_type == "Polygon":
        return [coordinates[0]]
    if geom_type == "MultiPolygon":
        return [polygon[0] for polygon in coordinates]
    return []


def build_map_patches(geojson: dict, values_by_ibge: dict[int, float]):
    from matplotlib.patches import Polygon

    patches: list[Polygon] = []
    values: list[float] = []
    missing_features: list[int] = []

    for feature in geojson["features"]:
        props = feature.get("properties", {})
        ibge_code = int(props["GEOCODIGO"])
        value = values_by_ibge.get(ibge_code)
        if value is None or not math.isfinite(value):
            missing_features.append(ibge_code)
            continue

        for ring in polygon_rings(feature["geometry"]):
            patches.append(Polygon(ring, closed=True))
            values.append(value)

    if not patches:
        raise ValueError("No mappable municipality polygons were matched to incidence values.")

    return patches, values, missing_features


def save_map(
    incidence: pd.DataFrame,
    geojson_file: Path,
    output_png: Path,
    output_pdf: Path,
    output_tiff: Path,
    show_title: bool,
    start_year: int,
    end_year: int,
) -> None:
    cache_dir = output_png.parent / ".plot-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    from matplotlib.colors import Normalize
    from matplotlib.ticker import FuncFormatter

    with geojson_file.open(encoding="utf-8") as f:
        geojson = json.load(f)

    values_by_ibge = dict(
        zip(
            incidence["ibge_code"].astype(int),
            incidence["mean_annual_incidence_per_100k"].astype(float),
            strict=True,
        )
    )
    patches, values, missing_features = build_map_patches(geojson, values_by_ibge)

    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 600,
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, ax = plt.subplots(figsize=(6.6, 4.9), constrained_layout=True)
    norm = Normalize(vmin=min(values), vmax=max(values))
    collection = PatchCollection(
        patches,
        cmap="YlOrRd",
        norm=norm,
        edgecolor="#ffffff",
        linewidth=0.22,
    )
    collection.set_array(values)
    ax.add_collection(collection)
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.axis("off")

    if show_title:
        ax.set_title(
            f"Mean Annual Dengue Incidence, {start_year}-{end_year}",
            loc="left",
            pad=8,
        )

    cbar = fig.colorbar(collection, ax=ax, fraction=0.032, pad=0.015)
    cbar.outline.set_visible(False)
    cbar.set_label("Mean annual incidence\nper 100,000 residents", fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))

    fig.savefig(output_png, bbox_inches="tight")
    fig.savefig(output_pdf, bbox_inches="tight")
    fig.savefig(output_tiff, bbox_inches="tight")
    plt.close(fig)

    if missing_features:
        missing_unique = sorted(set(missing_features))
        print(f"GeoJSON features without incidence values: {len(missing_unique)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Map municipal mean annual dengue incidence in Rio de Janeiro state."
    )
    parser.add_argument("--project-dir", default=None)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--start-year", type=int, default=DEFAULT_START_YEAR)
    parser.add_argument("--end-year", type=int, default=DEFAULT_END_YEAR)
    parser.add_argument("--show-title", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project_dir = find_project_dir(args.project_dir)
    data_dir = project_dir / "data"
    output_dir = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else project_dir / "outputs" / "descriptive_figures"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    incidence = load_mean_annual_incidence(data_dir, args.start_year, args.end_year)
    csv_path = output_dir / "municipal_mean_annual_dengue_incidence.csv"
    png_path = output_dir / "municipal_mean_annual_dengue_incidence.png"
    pdf_path = output_dir / "municipal_mean_annual_dengue_incidence.pdf"
    tiff_path = output_dir / "municipal_mean_annual_dengue_incidence.tiff"

    incidence.to_csv(csv_path, index=False)
    save_map(
        incidence,
        data_dir / "RJ.json",
        png_path,
        pdf_path,
        tiff_path,
        args.show_title,
        args.start_year,
        args.end_year,
    )

    print(f"Project directory: {project_dir}")
    print(f"Municipalities mapped: {len(incidence)}")
    print(f"Years summarized: {args.start_year}-{args.end_year}")
    print(
        "Mean annual incidence range: "
        f"{incidence['mean_annual_incidence_per_100k'].min():.1f}-"
        f"{incidence['mean_annual_incidence_per_100k'].max():.1f} per 100,000"
    )
    print(f"Wrote {csv_path}")
    print(f"Wrote {png_path}")
    print(f"Wrote {pdf_path}")
    print(f"Wrote {tiff_path}")


if __name__ == "__main__":
    main()
