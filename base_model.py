import os
import re
import unicodedata
import warnings

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymc as pm
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=FutureWarning)


# =========================================================
# Configuration
# =========================================================
DATA_START_YEAR = 2017
DATA_END_YEAR = 2023

MAKE_PLOTS = True
SAVE_OUTPUTS = True

DRAWS = 1000
TUNE = 2000
CHAINS = 4
CORES = 4
TARGET_ACCEPT = 0.98
RANDOM_SEED = 42

BASE_COVARIATES = [
    "rainfall",
    "humidity",
    "temperature",
    "idhm",
]


# =========================================================
# Paths
# =========================================================
BASE_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(BASE_DIR, "data")
COMBINED_FILE = os.path.join(DATA_DIR, "complete_combined_datasets.csv")


# =========================================================
# Helpers
# =========================================================
def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = df.columns.str.strip().str.lower()
    return df


def normalize_name(name):
    if pd.isna(name):
        return np.nan

    name = str(name).strip().lower()
    name = re.sub(r"/[a-z]{2}$", "", name)
    name = unicodedata.normalize("NFKD", name)
    name = "".join(c for c in name if not unicodedata.combining(c))
    name = re.sub(r"\s+", " ", name).strip()
    return name


def iso_week_to_date(year_series, week_series):
    return pd.to_datetime(
        year_series.astype(str)
        + "-W"
        + week_series.astype(int).astype(str).str.zfill(2)
        + "-1",
        format="%G-W%V-%u",
        errors="coerce",
    )


def compute_metrics(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)

    mae = float(np.mean(np.abs(y_true - y_pred)))
    rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
    wape = float(np.sum(np.abs(y_true - y_pred)) / max(np.sum(np.abs(y_true)), 1e-9))
    accuracy_pct = max(0.0, 100.0 * (1.0 - wape))

    sst = float(np.sum((y_true - y_true.mean()) ** 2))
    sse = float(np.sum((y_true - y_pred) ** 2))
    r2 = float(1.0 - sse / sst) if sst > 0 else np.nan

    return {
        "mae": mae,
        "rmse": rmse,
        "wape": wape,
        "accuracy_pct": accuracy_pct,
        "r2": r2,
    }


def da_mean(x):
    return float(np.asarray(x).mean())


def summarize_diagnostics(trace):
    rhat = az.rhat(trace)
    ess_bulk = az.ess(trace, method="bulk")
    ess_tail = az.ess(trace, method="tail")

    variables = ["beta", "alpha_municipio", "sigma_municipio"]
    return {
        "mean_rhat": float(np.mean([da_mean(rhat[v]) for v in variables])),
        "mean_ess_bulk": float(np.mean([da_mean(ess_bulk[v]) for v in variables])),
        "mean_ess_tail": float(np.mean([da_mean(ess_tail[v]) for v in variables])),
    }


# =========================================================
# Data preparation
# =========================================================
def build_model_dataframe() -> pd.DataFrame:
    df = clean_columns(pd.read_csv(COMBINED_FILE))
    original_rows = len(df)

    df["municipio"] = df["municipio"].apply(normalize_name)

    numeric_cols = ["year", "week", "cases", *BASE_COVARIATES]
    if "population" in df.columns:
        numeric_cols.append("population")

    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    core_required = ["municipio", "year", "week", "cases"]
    rows_before_core_drop = len(df)
    df = df.dropna(subset=core_required).copy()
    core_drop_count = rows_before_core_drop - len(df)

    df["year"] = df["year"].astype(int)
    df["week"] = df["week"].astype(int)

    rows_before_year_filter = len(df)
    df = df[(df["year"] >= DATA_START_YEAR) & (df["year"] <= DATA_END_YEAR)].copy()
    year_filter_drop_count = rows_before_year_filter - len(df)

    df["date"] = iso_week_to_date(df["year"], df["week"])
    rows_before_date_drop = len(df)
    df = df.dropna(subset=["date"]).copy()
    date_drop_count = rows_before_date_drop - len(df)

    rows_before_covariate_drop = len(df)
    df = df.dropna(subset=BASE_COVARIATES).copy()
    covariate_drop_count = rows_before_covariate_drop - len(df)

    df["cases"] = df["cases"].clip(lower=0).astype(int)
    df = df.sort_values(["municipio", "date"]).reset_index(drop=True)

    print("Missing-data handling: rows with NaN in required fields are dropped.")
    print(f"Original rows: {original_rows}")
    print(f"Dropped missing municipio/year/week/cases: {core_drop_count}")
    print(f"Dropped outside {DATA_START_YEAR}-{DATA_END_YEAR}: {year_filter_drop_count}")
    print(f"Dropped invalid ISO week dates: {date_drop_count}")
    print(f"Dropped missing covariates {BASE_COVARIATES}: {covariate_drop_count}")
    print("Model rows:", len(df))
    print("Municipios:", df["municipio"].nunique())
    print("Years:", sorted(df["year"].unique().tolist()))
    print("Covariates:", BASE_COVARIATES)

    return df


# =========================================================
# Model
# =========================================================
def prepare_arrays(df: pd.DataFrame):
    df = df.copy()
    df["municipio_idx"], municipios = pd.factorize(df["municipio"], sort=True)

    week_levels = sorted(df["week"].unique().tolist())
    week_to_idx = {week: idx for idx, week in enumerate(week_levels)}
    df["week_idx"] = df["week"].map(week_to_idx).astype(int)

    year_levels = sorted(df["year"].unique().tolist())
    year_to_idx = {year: idx for idx, year in enumerate(year_levels)}
    df["year_idx"] = df["year"].map(year_to_idx).astype(int)

    scaler = StandardScaler()
    X = scaler.fit_transform(df[BASE_COVARIATES].to_numpy(dtype=float))

    return {
        "df": df,
        "X": X,
        "y": df["cases"].to_numpy(dtype=int),
        "municipio_idx": df["municipio_idx"].to_numpy(dtype=int),
        "week_idx": df["week_idx"].to_numpy(dtype=int),
        "year_idx": df["year_idx"].to_numpy(dtype=int),
        "municipios": municipios,
        "week_levels": week_levels,
        "year_levels": year_levels,
        "scaler": scaler,
    }


def fit_model(inputs):
    coords = {
        "obs_id": np.arange(len(inputs["y"])),
        "covariate": BASE_COVARIATES,
        "municipio": inputs["municipios"],
        "week": inputs["week_levels"],
        "year": inputs["year_levels"],
    }

    with pm.Model(coords=coords) as model:
        X = pm.Data("X", inputs["X"], dims=("obs_id", "covariate"))
        municipio_idx = pm.Data("municipio_idx", inputs["municipio_idx"], dims="obs_id")
        week_idx = pm.Data("week_idx", inputs["week_idx"], dims="obs_id")
        year_idx = pm.Data("year_idx", inputs["year_idx"], dims="obs_id")

        intercept = pm.Normal("intercept", mu=0, sigma=2)
        beta = pm.Normal("beta", mu=0, sigma=1, dims="covariate")

        sigma_municipio = pm.Exponential("sigma_municipio", 1.0)
        alpha_municipio_raw = pm.Normal(
            "alpha_municipio_raw",
            mu=0,
            sigma=1,
            dims="municipio",
        )
        alpha_municipio = pm.Deterministic(
            "alpha_municipio",
            alpha_municipio_raw * sigma_municipio,
            dims="municipio",
        )

        sigma_week = pm.Exponential("sigma_week", 1.0)
        week_effect = pm.Normal("week_effect", mu=0, sigma=sigma_week, dims="week")

        sigma_year = pm.Exponential("sigma_year", 1.0)
        year_effect = pm.Normal("year_effect", mu=0, sigma=sigma_year, dims="year")

        eta = (
            intercept
            + alpha_municipio[municipio_idx]
            + week_effect[week_idx]
            + year_effect[year_idx]
            + pm.math.dot(X, beta)
        )
        mu = pm.Deterministic("mu", pm.math.exp(eta), dims="obs_id")

        alpha_nb = pm.Exponential("alpha_nb", 1.0)
        pm.NegativeBinomial(
            "cases",
            mu=mu,
            alpha=alpha_nb,
            observed=inputs["y"],
            dims="obs_id",
        )

        trace = pm.sample(
            draws=DRAWS,
            tune=TUNE,
            chains=CHAINS,
            cores=CORES,
            target_accept=TARGET_ACCEPT,
            random_seed=RANDOM_SEED,
            return_inferencedata=True,
        )
        posterior_predictive = pm.sample_posterior_predictive(
            trace,
            var_names=["cases"],
            random_seed=RANDOM_SEED,
            return_inferencedata=False,
        )

    return model, trace, posterior_predictive


def main():
    df = build_model_dataframe()
    inputs = prepare_arrays(df)
    _, trace, posterior_predictive = fit_model(inputs)

    pred_mean = posterior_predictive["cases"].mean(axis=(0, 1))
    metrics = compute_metrics(inputs["y"], pred_mean)
    diagnostics = summarize_diagnostics(trace)

    print("\nModel metrics:")
    for key, value in metrics.items():
        print(f"{key}: {value:.4f}")

    print("\nDiagnostics:")
    for key, value in diagnostics.items():
        print(f"{key}: {value:.4f}")

    print("\nPosterior summary:")
    print(
        az.summary(
            trace,
            var_names=[
                "intercept",
                "beta",
                "sigma_municipio",
                "sigma_week",
                "sigma_year",
                "alpha_nb",
            ],
        )
    )

    if SAVE_OUTPUTS:
        metrics_path = os.path.join(BASE_DIR, "base_model_metrics.csv")
        predictions_path = os.path.join(BASE_DIR, "base_model_predictions.csv")

        pd.DataFrame([{**metrics, **diagnostics}]).to_csv(metrics_path, index=False)

        output_df = inputs["df"][["municipio", "year", "week", "date", "cases"]].copy()
        output_df["predicted_cases"] = pred_mean
        output_df.to_csv(predictions_path, index=False)

        print(f"\nSaved metrics to: {metrics_path}")
        print(f"Saved predictions to: {predictions_path}")

    if MAKE_PLOTS:
        plt.figure(figsize=(8, 8))
        plt.scatter(inputs["y"], pred_mean, alpha=0.2)
        line_min = min(float(inputs["y"].min()), float(pred_mean.min()))
        line_max = max(float(inputs["y"].max()), float(pred_mean.max()))
        plt.plot([line_min, line_max], [line_min, line_max], "r--")
        plt.xlabel("Actual cases")
        plt.ylabel("Posterior predictive mean cases")
        plt.title("Base hierarchical model: actual vs predicted")
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    main()
