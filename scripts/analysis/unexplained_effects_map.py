"""Create the S2 unexplained spatial-effects map.

This Python script is a convenience wrapper around the S2 R-INLA spatial-effect
extraction. S2 is fit in R because it uses INLA/BYM2, then this script reads the
resulting CSV and maps residual spatial relative risk using data/RJ.json.

Run from the project root:

    python unexplained_effects_map.py

Outputs:
    outputs/s2_unexplained_spatial_effects.csv
    outputs/s2_unexplained_spatial_relative_risk_map.png
"""

from __future__ import annotations

import csv
import json
import math
import os
import subprocess
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path.cwd() / "outputs" / ".matplotlib-cache"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path.cwd() / "outputs" / ".plot-cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Polygon


def find_project_dir() -> Path:
    override = os.environ.get("HBM_PROJECT_DIR")
    if override:
        return Path(override).expanduser().resolve()

    here = Path(__file__).resolve().parent
    for candidate in [Path.cwd(), here, *here.parents]:
        if (candidate / "data" / "complete_combined_datasets.csv").exists():
            return candidate.resolve()

    raise FileNotFoundError(
        "Could not find project directory. Run from the project root or set "
        "HBM_PROJECT_DIR=/path/to/Hierarchical-Bayesian-Model."
    )


PROJECT_DIR = find_project_dir()
DATA_DIR = PROJECT_DIR / "data"
SPATIAL_R_DIR = PROJECT_DIR / "models" / "r_inla" / "spatial"
OUTPUT_DIR = PROJECT_DIR / "outputs"

RJ_GEOJSON = DATA_DIR / "RJ.json"
EFFECTS_CSV = OUTPUT_DIR / "s2_unexplained_spatial_effects.csv"
MAP_PNG = OUTPUT_DIR / "s2_unexplained_spatial_relative_risk_map.png"
S2_EFFECT_SCRIPT = SPATIAL_R_DIR / "map_s2_unexplained_effects.R"


def run_s2_effect_extraction() -> None:
    """Run the R-INLA S2 extraction script if the effects CSV is unavailable."""
    if not S2_EFFECT_SCRIPT.exists():
        raise FileNotFoundError(
            f"Missing {S2_EFFECT_SCRIPT}. Run/create the R S2 effects script first."
        )

    subprocess.run(
        ["Rscript", str(S2_EFFECT_SCRIPT)],
        cwd=PROJECT_DIR,
        check=True,
    )


def load_effects() -> dict[int, dict[str, str]]:
    if not EFFECTS_CSV.exists():
        run_s2_effect_extraction()

    effects: dict[int, dict[str, str]] = {}
    with EFFECTS_CSV.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required = {"ibge_code", "residual_spatial_rr"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{EFFECTS_CSV} is missing columns: {sorted(missing)}")

        for row in reader:
            effects[int(row["ibge_code"])] = row

    return effects


def polygon_rings(geometry: dict) -> list[list[list[float]]]:
    """Return exterior rings from Polygon or MultiPolygon GeoJSON geometry."""
    geom_type = geometry.get("type")
    coordinates = geometry.get("coordinates", [])

    if geom_type == "Polygon":
        return [coordinates[0]]
    if geom_type == "MultiPolygon":
        return [polygon[0] for polygon in coordinates]

    return []


def build_map_patches(geojson: dict, effects: dict[int, dict[str, str]]):
    patches: list[Polygon] = []
    values: list[float] = []
    missing_features: list[int] = []

    for feature in geojson["features"]:
        props = feature.get("properties", {})
        ibge_code = int(props["GEOCODIGO"])
        effect = effects.get(ibge_code)

        if effect is None:
            missing_features.append(ibge_code)
            continue

        value = float(effect["residual_spatial_rr"])
        if not math.isfinite(value):
            missing_features.append(ibge_code)
            continue

        for ring in polygon_rings(feature["geometry"]):
            patches.append(Polygon(ring, closed=True))
            values.append(value)

    if not patches:
        raise ValueError("No mappable municipality polygons were matched to S2 effects.")

    return patches, values, missing_features


def relative_risk_norm(values: list[float]):
    vmin = min(values)
    vmax = max(values)

    if vmin < 1.0 < vmax:
        return TwoSlopeNorm(vmin=vmin, vcenter=1.0, vmax=vmax)

    if vmin == vmax:
        padding = max(abs(vmin) * 0.05, 0.05)
        return Normalize(vmin=vmin - padding, vmax=vmax + padding)

    return Normalize(vmin=vmin, vmax=vmax)


def plot_unexplained_effects() -> None:
    effects = load_effects()

    with RJ_GEOJSON.open(encoding="utf-8") as f:
        geojson = json.load(f)

    patches, values, missing_features = build_map_patches(geojson, effects)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(12, 8))
    fig.subplots_adjust(top=0.84, bottom=0.12)
    norm = relative_risk_norm(values)
    collection = PatchCollection(
        patches,
        cmap="coolwarm",
        norm=norm,
        edgecolor="#bdbdbd",
        linewidth=0.35,
    )
    collection.set_array(values)
    ax.add_collection(collection)
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.axis("off")

    cbar = fig.colorbar(collection, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_ticks([norm.vmin, 1.0, norm.vmax])
    cbar.set_ticklabels([f"{norm.vmin:.2f}\nLower", "1.0\nNo change", f"{norm.vmax:.2f}\nHigher"])
    cbar.set_label("Residual spatial RR")

    fig.text(
        0.03,
        0.955,
        "S2 Residual Spatial Relative Risk",
        fontsize=18,
        weight="bold",
        ha="left",
        va="top",
    )
    fig.text(
        0.03,
        0.915,
        "exp(posterior mean BYM2 spatial effect); adjusted for lagged weather, IDHM, own-case lag, and neighboring-case lag",
        fontsize=11,
        ha="left",
        va="top",
    )
    fig.text(
        0.03,
        0.055,
        "Blue values below 1 indicate lower unexplained spatial dengue risk; white at RR = 1 indicates no change; red values above 1 indicate higher unexplained spatial risk.",
        fontsize=10,
        ha="left",
        va="bottom",
    )

    fig.tight_layout()
    fig.savefig(MAP_PNG, dpi=300, bbox_inches="tight")

    print(f"Wrote map: {MAP_PNG}")
    print(f"Used effects: {EFFECTS_CSV}")
    if missing_features:
        print(f"GeoJSON features without mapped S2 effects: {len(missing_features)}")


if __name__ == "__main__":
    plot_unexplained_effects()
