#!/usr/bin/env python3
"""
Create a publication-ready map of Rio de Janeiro IBGE intermediate regions.

The map uses data/RJ.json for municipality boundaries and joins those polygons
to data/rj_ibge_intermediate_regions.csv by IBGE municipality code.

Default outputs:
  outputs/descriptive_figures/rj_intermediate_regions_map.png
  outputs/descriptive_figures/rj_intermediate_regions_map.pdf
  outputs/descriptive_figures/rj_intermediate_regions_map.tiff
  outputs/descriptive_figures/rj_geographic_regions_reference_style_map.png
  outputs/descriptive_figures/rj_geographic_regions_reference_style_map.pdf
  outputs/descriptive_figures/rj_geographic_regions_reference_style_map.tiff
  outputs/descriptive_figures/rj_intermediate_regions_journal_map.png
  outputs/descriptive_figures/rj_intermediate_regions_journal_map.pdf
  outputs/descriptive_figures/rj_intermediate_regions_journal_map.tiff

Run from the project root:
  python scripts/figures/plot_rj_intermediate_regions_map.py
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import pandas as pd


REGION_FIELD = "regiao_intermediaria_nome"
OUTPUT_STEM = "rj_intermediate_regions_map"
REFERENCE_STYLE_OUTPUT_STEM = "rj_geographic_regions_reference_style_map"
JOURNAL_OUTPUT_STEM = "rj_intermediate_regions_journal_map"

REGION_COLORS = {
    "Rio de Janeiro": "#4f7cac",
    "Volta Redonda - Barra Mansa": "#c27d38",
    "Petrópolis": "#5f8f7a",
    "Campos dos Goytacazes": "#9a4f5c",
    "Macaé - Rio das Ostras - Cabo Frio": "#7b6fa6",
}

BASE_REGION_COLORS = {
    "Rio de Janeiro": "#c99466",
    "Volta Redonda - Barra Mansa": "#b77a3d",
    "Petrópolis": "#a66d52",
    "Campos dos Goytacazes": "#c7773f",
    "Macaé - Rio das Ostras - Cabo Frio": "#8f796d",
}

JOURNAL_REGION_COLORS = {
    "Rio de Janeiro": "#8aa3b8",
    "Volta Redonda - Barra Mansa": "#c59a8f",
    "Petrópolis": "#8fa98d",
    "Campos dos Goytacazes": "#b493a6",
    "Macaé - Rio das Ostras - Cabo Frio": "#c8bd87",
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
        if (candidate / "data" / "RJ.json").exists():
            return candidate.resolve()

    raise FileNotFoundError(
        "Could not find project directory. Run from the project root or set "
        "HBM_PROJECT_DIR=/path/to/Hierarchical-Bayesian-Model."
    )


def polygon_rings(geometry: dict) -> list[list[list[float]]]:
    geom_type = geometry.get("type")
    coordinates = geometry.get("coordinates", [])
    if geom_type == "Polygon":
        return [coordinates[0]]
    if geom_type == "MultiPolygon":
        return [polygon[0] for polygon in coordinates]
    return []


def load_region_lookup(data_dir: Path) -> pd.DataFrame:
    path = data_dir / "rj_ibge_intermediate_regions.csv"
    regions = pd.read_csv(path)
    required = {"ibge_code", "municipio", REGION_FIELD}
    missing = sorted(required.difference(regions.columns))
    if missing:
        raise ValueError(f"Missing required columns in {path}: {missing}")

    regions["ibge_code"] = pd.to_numeric(regions["ibge_code"], errors="coerce").astype("Int64")
    regions = regions.dropna(subset=["ibge_code", REGION_FIELD]).copy()
    regions["ibge_code"] = regions["ibge_code"].astype(int)
    keep = [
        "ibge_code",
        "municipio",
        "regiao_imediata_id",
        "regiao_imediata_nome",
        "regiao_intermediaria_id",
        REGION_FIELD,
    ]
    return regions[keep].drop_duplicates(subset=["ibge_code"])


def blend_with_white(hex_color: str, amount: float) -> str:
    hex_color = hex_color.lstrip("#")
    red = int(hex_color[0:2], 16)
    green = int(hex_color[2:4], 16)
    blue = int(hex_color[4:6], 16)
    blended = [round(channel + (255 - channel) * amount) for channel in (red, green, blue)]
    return "#{:02x}{:02x}{:02x}".format(*blended)


def format_degree(value: float, suffix: str) -> str:
    return f"{abs(value):.0f} deg {suffix}"


def build_region_patches(geojson: dict, regions: pd.DataFrame):
    from matplotlib.patches import Polygon

    region_by_ibge = dict(
        zip(regions["ibge_code"].astype(int), regions[REGION_FIELD].astype(str), strict=True)
    )
    patches: list[Polygon] = []
    colors: list[str] = []
    missing_features: list[int] = []

    for feature in geojson["features"]:
        props = feature.get("properties", {})
        ibge_code = int(props["GEOCODIGO"])
        region = region_by_ibge.get(ibge_code)
        if region is None:
            missing_features.append(ibge_code)
            continue

        color = REGION_COLORS.get(region, "#b0b0b0")
        for ring in polygon_rings(feature["geometry"]):
            patches.append(Polygon(ring, closed=True))
            colors.append(color)

    if not patches:
        raise ValueError("No mappable municipality polygons were matched to region values.")

    return patches, colors, missing_features


def add_region_labels(ax, geojson: dict, regions: pd.DataFrame) -> None:
    region_by_ibge = dict(
        zip(regions["ibge_code"].astype(int), regions[REGION_FIELD].astype(str), strict=True)
    )
    coords_by_region: dict[str, list[tuple[float, float]]] = {}
    for feature in geojson["features"]:
        ibge_code = int(feature.get("properties", {}).get("GEOCODIGO"))
        region = region_by_ibge.get(ibge_code)
        if region is None:
            continue
        for ring in polygon_rings(feature["geometry"]):
            coords_by_region.setdefault(region, []).extend((float(x), float(y)) for x, y in ring)

    label_offsets = {
        "Rio de Janeiro": (0.12, -0.12),
        "Volta Redonda - Barra Mansa": (-0.04, -0.10),
        "Petrópolis": (0.0, 0.06),
        "Campos dos Goytacazes": (0.10, 0.02),
        "Macaé - Rio das Ostras - Cabo Frio": (0.16, -0.02),
    }
    label_text = {
        "Rio de Janeiro": "Rio de\nJaneiro",
        "Volta Redonda - Barra Mansa": "Volta Redonda\nBarra Mansa",
        "Petrópolis": "Petrópolis",
        "Campos dos Goytacazes": "Campos dos\nGoytacazes",
        "Macaé - Rio das Ostras - Cabo Frio": "Macaé - Rio das Ostras\nCabo Frio",
    }

    for region, coords in coords_by_region.items():
        xs = [point[0] for point in coords if math.isfinite(point[0])]
        ys = [point[1] for point in coords if math.isfinite(point[1])]
        if not xs or not ys:
            continue
        x = (min(xs) + max(xs)) / 2
        y = (min(ys) + max(ys)) / 2
        dx, dy = label_offsets.get(region, (0, 0))
        ax.text(
            x + dx,
            y + dy,
            label_text.get(region, region),
            ha="center",
            va="center",
            fontsize=7.5,
            color="#222222",
            linespacing=1.05,
        )


def save_region_map(project_dir: Path, output_dir: Path, show_title: bool) -> None:
    cache_dir = output_dir / ".plot-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Patch

    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 600,
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 12,
            "axes.labelsize": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    with (project_dir / "data" / "RJ.json").open(encoding="utf-8") as f:
        geojson = json.load(f)
    regions = load_region_lookup(project_dir / "data")
    patches, colors, missing_features = build_region_patches(geojson, regions)

    fig, ax = plt.subplots(figsize=(8.0, 5.7))
    fig.patch.set_facecolor("white")
    collection = PatchCollection(
        patches,
        facecolor=colors,
        edgecolor="white",
        linewidth=0.45,
        antialiased=True,
    )
    ax.add_collection(collection)

    outline = PatchCollection(
        patches,
        facecolor="none",
        edgecolor="#3f3f3f",
        linewidth=0.18,
        antialiased=True,
    )
    ax.add_collection(outline)
    add_region_labels(ax, geojson, regions)

    all_x = [xy[0] for patch in patches for xy in patch.get_xy()]
    all_y = [xy[1] for patch in patches for xy in patch.get_xy()]
    ax.set_xlim(min(all_x) - 0.12, max(all_x) + 0.12)
    ax.set_ylim(min(all_y) - 0.12, max(all_y) + 0.12)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")

    if show_title:
        ax.set_title("Rio de Janeiro IBGE Intermediate Regions", loc="left", pad=8)

    legend_order = [
        "Rio de Janeiro",
        "Volta Redonda - Barra Mansa",
        "Petrópolis",
        "Campos dos Goytacazes",
        "Macaé - Rio das Ostras - Cabo Frio",
    ]
    handles = [
        Patch(facecolor=REGION_COLORS[region], edgecolor="white", label=region)
        for region in legend_order
    ]
    legend = ax.legend(
        handles=handles,
        title="IBGE intermediate region",
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.08),
        borderaxespad=0,
        handlelength=1.2,
        handleheight=1.0,
        labelspacing=0.45,
        columnspacing=1.35,
        ncol=2,
        fontsize=8.0,
        title_fontsize=8.8,
    )
    legend._legend_box.align = "left"

    fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.17)
    output_stem = output_dir / OUTPUT_STEM
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_stem.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".tiff"), bbox_inches="tight")
    plt.close(fig)

    if missing_features:
        missing = ", ".join(str(code) for code in sorted(missing_features))
        print(f"Warning: GeoJSON features missing region lookup: {missing}")
    print(f"Mapped municipalities: {regions['ibge_code'].nunique()}")
    print(f"Mapped regions: {regions[REGION_FIELD].nunique()}")
    print(f"Wrote map to: {output_stem.with_suffix('.png')}")


def add_reference_labels(ax, immediate_gdf, intermediate_gdf) -> None:
    immediate_labels = {
        "Rio de Janeiro": "Rio de\nJaneiro",
        "Angra dos Reis": "Angra\ndos Reis",
        "Rio Bonito": "Rio\nBonito",
        "Volta Redonda - Barra Mansa": "Volta Redonda\nBarra Mansa",
        "Resende": "Resende",
        "Valença": "Valenca",
        "Petrópolis": "Petropolis",
        "Nova Friburgo": "Nova\nFriburgo",
        "Três Rios - Paraíba do Sul": "Tres Rios\nParaiba do Sul",
        "Campos dos Goytacazes": "Campos dos\nGoytacazes",
        "Itaperuna": "Itaperuna",
        "Santo Antônio de Pádua": "Santo Antonio\nde Padua",
        "Cabo Frio": "Cabo Frio",
        "Macaé - Rio das Ostras": "Macae\nRio das Ostras",
    }
    immediate_offsets = {
        "Rio de Janeiro": (0.02, -0.03),
        "Angra dos Reis": (-0.06, -0.06),
        "Rio Bonito": (0.0, -0.04),
        "Volta Redonda - Barra Mansa": (-0.03, 0.0),
        "Resende": (-0.04, 0.0),
        "Valença": (0.02, 0.0),
        "Petrópolis": (0.0, -0.02),
        "Nova Friburgo": (0.06, 0.02),
        "Três Rios - Paraíba do Sul": (-0.02, 0.02),
        "Campos dos Goytacazes": (0.08, 0.0),
        "Itaperuna": (-0.04, 0.04),
        "Santo Antônio de Pádua": (-0.05, -0.01),
        "Cabo Frio": (0.05, -0.01),
        "Macaé - Rio das Ostras": (0.02, 0.01),
    }
    for _, row in immediate_gdf.iterrows():
        point = row.geometry.representative_point()
        name = row["regiao_imediata_nome"]
        dx, dy = immediate_offsets.get(name, (0, 0))
        ax.text(
            point.x + dx,
            point.y + dy,
            immediate_labels.get(name, name),
            ha="center",
            va="center",
            fontsize=4.6,
            color="#342922",
            alpha=0.82,
            linespacing=0.95,
            zorder=7,
        )

    intermediate_labels = {
        "Rio de Janeiro": "Rio de Janeiro",
        "Volta Redonda - Barra Mansa": "Volta Redonda - Barra Mansa",
        "Petrópolis": "Petropolis",
        "Campos dos Goytacazes": "Campos dos Goytacazes",
        "Macaé - Rio das Ostras - Cabo Frio": "Macae - Rio das Ostras\nCabo Frio",
    }
    intermediate_offsets = {
        "Rio de Janeiro": (0.05, 0.02),
        "Volta Redonda - Barra Mansa": (-0.10, -0.02),
        "Petrópolis": (-0.02, 0.15),
        "Campos dos Goytacazes": (0.04, 0.20),
        "Macaé - Rio das Ostras - Cabo Frio": (0.18, -0.08),
    }
    for _, row in intermediate_gdf.iterrows():
        point = row.geometry.representative_point()
        name = row[REGION_FIELD]
        dx, dy = intermediate_offsets.get(name, (0, 0))
        ax.text(
            point.x + dx,
            point.y + dy,
            intermediate_labels.get(name, name),
            ha="center",
            va="center",
            fontsize=8.2,
            color="#df2f2f",
            fontweight="bold",
            zorder=8,
        )


def add_text_box(ax, x: float, y: float, title: str, lines: list[str], fontsize: float = 5.8) -> None:
    text = title + "\n" + "\n".join(lines)
    ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=fontsize,
        color="#222222",
        linespacing=1.15,
        bbox={
            "boxstyle": "square,pad=0.35",
            "facecolor": "white",
            "edgecolor": "#777777",
            "linewidth": 0.7,
            "alpha": 0.92,
        },
        zorder=20,
    )


def add_reference_style_map(project_dir: Path, output_dir: Path, show_title: bool) -> None:
    cache_dir = output_dir / ".plot-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))

    import geopandas as gpd
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 600,
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 12,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    regions = load_region_lookup(project_dir / "data")
    gdf = gpd.read_file(project_dir / "data" / "RJ.json")
    gdf["ibge_code"] = pd.to_numeric(gdf["GEOCODIGO"], errors="coerce").astype("Int64")
    gdf = gdf.merge(regions, on="ibge_code", how="left")
    missing = gdf.loc[gdf[REGION_FIELD].isna(), "GEOCODIGO"].astype(str).tolist()
    if missing:
        raise ValueError(f"GeoJSON municipalities missing region lookup: {', '.join(missing)}")

    immediate_gdf = gdf.dissolve(
        by=["regiao_imediata_id", "regiao_imediata_nome", "regiao_intermediaria_id", REGION_FIELD],
        as_index=False,
    )
    intermediate_gdf = gdf.dissolve(by=REGION_FIELD, as_index=False)

    immediate_colors: dict[str, str] = {}
    for region_name, group in immediate_gdf.groupby(REGION_FIELD, sort=False):
        names = group.sort_values("regiao_imediata_id")["regiao_imediata_nome"].tolist()
        base = BASE_REGION_COLORS.get(region_name, "#b9895a")
        amounts = [0.06, 0.20, 0.34, 0.46]
        for index, name in enumerate(names):
            immediate_colors[name] = blend_with_white(base, amounts[index % len(amounts)])
    gdf["map_color"] = gdf["regiao_imediata_nome"].map(immediate_colors)

    minx, miny, maxx, maxy = gdf.total_bounds
    fig, ax = plt.subplots(figsize=(12.0, 8.5))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("#b9d8ef")

    land_x = [minx - 0.55, maxx + 0.55, maxx + 0.55, minx - 0.55]
    land_y = [miny + 0.30, miny + 0.30, maxy + 0.45, maxy + 0.45]
    ax.fill(land_x, land_y, color="#eeeeeb", zorder=0)

    gdf.plot(
        ax=ax,
        color=gdf["map_color"],
        edgecolor="#5a4535",
        linewidth=0.28,
        alpha=0.78,
        zorder=3,
    )
    immediate_gdf.boundary.plot(ax=ax, color="#3f3128", linewidth=0.45, alpha=0.8, zorder=5)
    intermediate_gdf.boundary.plot(ax=ax, color="#e02f2f", linewidth=1.1, zorder=6)
    gdf.dissolve().boundary.plot(ax=ax, color="#ff2d2d", linewidth=1.35, zorder=7)
    add_reference_labels(ax, immediate_gdf, intermediate_gdf)

    ax.set_xlim(minx - 0.35, maxx + 0.35)
    ax.set_ylim(miny - 0.75, maxy + 0.35)
    ax.set_aspect("equal", adjustable="box")
    xticks = [-45, -44, -43, -42, -41]
    yticks = [-23, -22, -21]
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels([format_degree(value, "W") for value in xticks], fontsize=5.8)
    ax.set_yticklabels([format_degree(value, "S") for value in yticks], fontsize=5.8)
    ax.tick_params(top=True, labeltop=True, right=True, labelright=True, length=2.5, colors="#444444")
    ax.grid(color="#7e9cb4", linewidth=0.45, alpha=0.45, zorder=1)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_edgecolor("#4d4d4d")

    ax.text(
        0.985,
        0.965,
        "Regioes Geograficas\nEstado do Rio de Janeiro",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=13,
        fontweight="bold",
        color="#222222",
        linespacing=1.05,
        bbox={
            "boxstyle": "square,pad=0.35",
            "facecolor": "white",
            "edgecolor": "#8a8a8a",
            "linewidth": 0.7,
            "alpha": 0.92,
        },
        zorder=22,
    )
    ax.text(-42.9, miny - 0.47, "O C E A N O   A T L A N T I C O", color="#2f79bd", fontsize=8, alpha=0.75)

    immediate_rows = (
        regions[["regiao_imediata_id", "regiao_imediata_nome"]]
        .drop_duplicates()
        .sort_values("regiao_imediata_id")
    )
    immediate_lines = [
        f"{int(row.regiao_imediata_id)} - {row.regiao_imediata_nome}"
        for row in immediate_rows.itertuples(index=False)
    ]
    add_text_box(ax, 0.018, 0.026, "Regioes Geograficas Imediatas", immediate_lines, fontsize=5.0)

    intermediate_rows = (
        regions[["regiao_intermediaria_id", REGION_FIELD]]
        .drop_duplicates()
        .sort_values("regiao_intermediaria_id")
    )
    intermediate_lines = [
        f"{int(row.regiao_intermediaria_id)} - {getattr(row, REGION_FIELD)}"
        for row in intermediate_rows.itertuples(index=False)
    ]
    add_text_box(ax, 0.305, 0.040, "Regioes Geograficas\nIntermediarias", intermediate_lines, fontsize=5.3)

    convention_handles = [
        Line2D([0], [0], color="#ff2d2d", linewidth=1.4, label="State/intermediate boundary"),
        Line2D([0], [0], color="#3f3128", linewidth=0.55, label="Immediate-region boundary"),
        Line2D([0], [0], color="#5a4535", linewidth=0.35, label="Municipality boundary"),
        Patch(facecolor="#b9d8ef", edgecolor="#4e8fc4", label="Water"),
    ]
    convention = ax.legend(
        handles=convention_handles,
        title="Conventions",
        loc="lower right",
        bbox_to_anchor=(0.76, 0.052),
        frameon=True,
        fancybox=False,
        framealpha=0.92,
        edgecolor="#777777",
        facecolor="white",
        fontsize=5.3,
        title_fontsize=5.8,
        handlelength=1.8,
        borderpad=0.45,
        labelspacing=0.35,
    )
    ax.add_artist(convention)

    scalebar_y = miny - 0.54
    scalebar_x = maxx - 1.25
    scalebar_width = 0.48
    ax.plot([scalebar_x, scalebar_x + scalebar_width], [scalebar_y, scalebar_y], color="#222222", linewidth=1.2)
    ax.plot([scalebar_x, scalebar_x], [scalebar_y - 0.025, scalebar_y + 0.025], color="#222222", linewidth=1.0)
    ax.plot(
        [scalebar_x + scalebar_width, scalebar_x + scalebar_width],
        [scalebar_y - 0.025, scalebar_y + 0.025],
        color="#222222",
        linewidth=1.0,
    )
    ax.text(scalebar_x + scalebar_width / 2, scalebar_y - 0.07, "50 km", ha="center", va="top", fontsize=5.8)

    if show_title:
        ax.set_title("Rio de Janeiro geographic regions", loc="left", pad=8)

    output_stem = output_dir / REFERENCE_STYLE_OUTPUT_STEM
    fig.savefig(output_stem.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".tiff"), bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote reference-style map to: {output_stem.with_suffix('.png')}")


def add_journal_style_map(project_dir: Path, output_dir: Path, show_title: bool) -> None:
    cache_dir = output_dir / ".plot-cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))

    import geopandas as gpd
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 600,
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 10,
            "axes.labelsize": 8,
            "legend.fontsize": 6.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    regions = load_region_lookup(project_dir / "data")
    gdf = gpd.read_file(project_dir / "data" / "RJ.json")
    gdf["ibge_code"] = pd.to_numeric(gdf["GEOCODIGO"], errors="coerce").astype("Int64")
    gdf = gdf.merge(regions, on="ibge_code", how="left")
    missing = gdf.loc[gdf[REGION_FIELD].isna(), "GEOCODIGO"].astype(str).tolist()
    if missing:
        raise ValueError(f"GeoJSON municipalities missing region lookup: {', '.join(missing)}")

    intermediate_gdf = gdf.dissolve(by=REGION_FIELD, as_index=False)
    gdf["map_color"] = gdf[REGION_FIELD].map(JOURNAL_REGION_COLORS)

    minx, miny, maxx, maxy = gdf.total_bounds
    fig, ax = plt.subplots(figsize=(7.2, 4.95))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("#fbfbfa")

    water_top = miny + 0.20
    ax.axhspan(miny - 0.34, water_top, color="#eef5f8", zorder=0)
    ax.axhspan(water_top, maxy + 0.18, color="#fbfbfa", zorder=0)

    gdf.plot(
        ax=ax,
        color=gdf["map_color"],
        edgecolor="#f4f4f2",
        linewidth=0.22,
        alpha=0.92,
        zorder=3,
    )
    intermediate_gdf.boundary.plot(ax=ax, color="#3f3f3f", linewidth=0.55, zorder=5)
    gdf.dissolve().boundary.plot(ax=ax, color="#252525", linewidth=0.75, zorder=6)

    label_offsets = {
        "Rio de Janeiro": (0.08, 0.06),
        "Volta Redonda - Barra Mansa": (-0.04, 0.00),
        "Petrópolis": (-0.05, 0.03),
        "Campos dos Goytacazes": (0.04, 0.08),
        "Macaé - Rio das Ostras - Cabo Frio": (0.14, -0.05),
    }
    label_text = {
        "Rio de Janeiro": "Rio de Janeiro",
        "Volta Redonda - Barra Mansa": "Volta Redonda-\nBarra Mansa",
        "Petrópolis": "Petropolis",
        "Campos dos Goytacazes": "Campos dos\nGoytacazes",
        "Macaé - Rio das Ostras - Cabo Frio": "Macae-Rio das Ostras-\nCabo Frio",
    }
    for _, row in intermediate_gdf.iterrows():
        point = row.geometry.representative_point()
        region = row[REGION_FIELD]
        dx, dy = label_offsets.get(region, (0, 0))
        ax.text(
            point.x + dx,
            point.y + dy,
            label_text.get(region, region),
            ha="center",
            va="center",
            fontsize=6.0,
            color="#2f2f2f",
            linespacing=1.0,
            zorder=8,
        )

    ax.text(
        -42.75,
        miny - 0.26,
        "Atlantic Ocean",
        color="#5e8daa",
        fontsize=6.8,
        alpha=0.72,
        ha="center",
        va="center",
    )

    ax.set_xlim(minx - 0.15, maxx + 0.15)
    ax.set_ylim(miny - 0.32, maxy + 0.16)
    ax.set_aspect("equal", adjustable="box")

    xticks = [-45, -44, -43, -42, -41]
    yticks = [-23, -22, -21]
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels([format_degree(value, "W") for value in xticks], fontsize=5.8)
    ax.set_yticklabels([format_degree(value, "S") for value in yticks], fontsize=5.8)
    ax.tick_params(top=True, labeltop=True, bottom=True, labelbottom=True, right=True, labelright=True, length=2.0)
    ax.grid(color="#b9c8d1", linewidth=0.35, alpha=0.42, zorder=1)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.55)
        spine.set_edgecolor("#6b6b6b")

    legend_order = [
        "Rio de Janeiro",
        "Volta Redonda - Barra Mansa",
        "Petrópolis",
        "Campos dos Goytacazes",
        "Macaé - Rio das Ostras - Cabo Frio",
    ]
    region_handles = [
        Patch(facecolor=JOURNAL_REGION_COLORS[region], edgecolor="none", label=region)
        for region in legend_order
    ]
    line_handles = [
        Line2D([0], [0], color="#252525", linewidth=0.75, label="State boundary"),
        Line2D([0], [0], color="#3f3f3f", linewidth=0.55, label="Intermediate-region boundary"),
        Line2D([0], [0], color="#cfcfca", linewidth=0.55, label="Municipality boundary"),
    ]
    legend = ax.legend(
        handles=[*region_handles, *line_handles],
        title="IBGE intermediate region",
        loc="upper left",
        bbox_to_anchor=(0.02, 0.975),
        ncol=1,
        frameon=True,
        fancybox=False,
        framealpha=0.90,
        edgecolor="#b0b0b0",
        facecolor="white",
        handlelength=1.35,
        labelspacing=0.35,
        borderpad=0.38,
        title_fontsize=7.2,
    )
    legend._legend_box.align = "left"

    scalebar_x = maxx - 1.08
    scalebar_y = miny - 0.20
    scalebar_width = 0.48
    ax.plot([scalebar_x, scalebar_x + scalebar_width], [scalebar_y, scalebar_y], color="#222222", linewidth=0.8)
    ax.plot([scalebar_x, scalebar_x], [scalebar_y - 0.016, scalebar_y + 0.016], color="#222222", linewidth=0.75)
    ax.plot(
        [scalebar_x + scalebar_width, scalebar_x + scalebar_width],
        [scalebar_y - 0.016, scalebar_y + 0.016],
        color="#222222",
        linewidth=0.75,
    )
    ax.text(scalebar_x + scalebar_width / 2, scalebar_y - 0.050, "50 km", ha="center", va="top", fontsize=5.7)

    ax.annotate(
        "N",
        xy=(0.965, 0.14),
        xytext=(0.965, 0.055),
        xycoords="axes fraction",
        textcoords="axes fraction",
        arrowprops={"arrowstyle": "-|>", "linewidth": 0.75, "color": "#222222"},
        ha="center",
        va="center",
        fontsize=6.5,
        color="#222222",
    )

    if show_title:
        ax.set_title("Rio de Janeiro IBGE intermediate regions", loc="left", pad=8)

    fig.subplots_adjust(left=0.05, right=0.98, top=0.94, bottom=0.08)
    output_stem = output_dir / JOURNAL_OUTPUT_STEM
    fig.savefig(output_stem.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".tiff"), bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote journal-style map to: {output_stem.with_suffix('.png')}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Map Rio de Janeiro IBGE intermediate regions.")
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
        else project_dir / "outputs" / "descriptive_figures"
    )
    save_region_map(project_dir, output_dir, args.show_title)
    add_reference_style_map(project_dir, output_dir, args.show_title)
    add_journal_style_map(project_dir, output_dir, args.show_title)


if __name__ == "__main__":
    main()
