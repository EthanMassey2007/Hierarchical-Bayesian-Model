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
TARGET_STATE = "RJ"
DATA_START_YEAR = 2017
DATA_END_YEAR = 2023

CASE_LAG_WEEKS = 2
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
    "log_population",
    "log_density_2022",
    "hub_degree",
    "hub_proximity_pct",
    "air_pass_in",
    "air_pass_out",
    "air_conec_in",
    "air_conec_out",
    "road_conec_in",
    "road_conec_out",
    "fluv_conec_in",
    "fluv_conec_out",
    "log_cases_lag",
]

SKEWED_FEATURES = [
    "air_pass_in",
    "air_pass_out",
    "air_conec_in",
    "air_conec_out",
    "road_conec_in",
    "road_conec_out",
    "fluv_conec_in",
    "fluv_conec_out",
]


# =========================================================
# Paths
# =========================================================
BASE_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(BASE_DIR, "data")

COMBINED_FILE = os.path.join(DATA_DIR, "complete_combined_datasets.csv")
MUNICIPIOS_FILE = os.path.join(DATA_DIR, "municipios.csv")
AERO_FILE = os.path.join(DATA_DIR, "aero_anac_2017_2023.parquet")
FLUVI_FILE = os.path.join(DATA_DIR, "fluvi_road_ibge.parquet")
HUB_FILE = os.path.join(DATA_DIR, "hub_pop_density.csv")


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
# Feature engineering
# =========================================================
def build_municipio_lookup() -> pd.DataFrame:
    municipios = clean_columns(pd.read_csv(MUNICIPIOS_FILE))
    municipios["state_code"] = municipios["city"].astype(str).str.extract(
        r"/([A-Z]{2})$",
        expand=False,
    )
    municipios = municipios[municipios["state_code"] == TARGET_STATE].copy()
    municipios["municipio"] = municipios["city"].apply(normalize_name)
    municipios["ibge_code"] = pd.to_numeric(municipios["ibgeid"], errors="coerce")

    municipios_lookup = (
        municipios[["municipio", "ibge_code"]]
        .dropna()
        .drop_duplicates()
        .assign(ibge_code=lambda d: d["ibge_code"].astype(int))
    )

    hub = clean_columns(pd.read_csv(HUB_FILE))
    hub = hub[hub["uf"].astype(str).str.upper() == TARGET_STATE].copy()
    hub["municipio"] = hub["nm_municipio"].apply(normalize_name)
    hub["ibge_code"] = pd.to_numeric(hub["co_ibge"], errors="coerce")
    hub_lookup = (
        hub[["municipio", "ibge_code"]]
        .dropna()
        .drop_duplicates()
        .assign(ibge_code=lambda d: d["ibge_code"].astype(int))
    )

    lookup = pd.concat([municipios_lookup, hub_lookup], ignore_index=True)
    duplicate_codes = lookup.groupby("municipio")["ibge_code"].nunique()
    ambiguous = duplicate_codes[duplicate_codes > 1]
    if not ambiguous.empty:
        raise ValueError(f"Ambiguous IBGE lookup entries: {ambiguous.to_dict()}")

    return lookup.drop_duplicates(subset=["municipio"], keep="first").reset_index(drop=True)


def build_hub_features() -> pd.DataFrame:
    hub = clean_columns(pd.read_csv(HUB_FILE))
    hub = hub[hub["uf"].astype(str).str.upper() == TARGET_STATE].copy()

    hub["municipio"] = hub["nm_municipio"].apply(normalize_name)
    hub["ibge_code"] = pd.to_numeric(hub["co_ibge"], errors="coerce")
    hub["hub_degree"] = pd.to_numeric(hub["grau"], errors="coerce")
    hub["hub_proximity_pct"] = pd.to_numeric(hub["ind_proxi_per"], errors="coerce")
    hub["density_2022"] = pd.to_numeric(hub["densidade_2022"], errors="coerce")
    hub["population_2022"] = pd.to_numeric(hub["populacao_2022"], errors="coerce")

    return (
        hub[
            [
                "municipio",
                "ibge_code",
                "hub_degree",
                "hub_proximity_pct",
                "density_2022",
                "population_2022",
            ]
        ]
        .dropna(subset=["municipio", "ibge_code"])
        .drop_duplicates(subset=["ibge_code"])
        .assign(ibge_code=lambda d: d["ibge_code"].astype(int))
    )


def build_aero_features(valid_ibge) -> pd.DataFrame:
    aero = clean_columns(pd.read_parquet(AERO_FILE))
    aero = aero.rename(columns={"ano": "year", "mes": "month"})

    numeric_cols = [
        "year",
        "month",
        "co_muni_ori",
        "co_muni_des",
        "aero_pass",
        "aero_pass_week",
        "aero_conec",
    ]
    for col in numeric_cols:
        aero[col] = pd.to_numeric(aero[col], errors="coerce")

    aero = aero.dropna(subset=["year", "month", "co_muni_ori", "co_muni_des"]).copy()
    for col in ["year", "month", "co_muni_ori", "co_muni_des"]:
        aero[col] = aero[col].astype(int)

    aero = aero[
        (aero["year"] >= DATA_START_YEAR)
        & (aero["year"] <= DATA_END_YEAR)
        & (
            aero["co_muni_ori"].isin(valid_ibge)
            | aero["co_muni_des"].isin(valid_ibge)
        )
    ].copy()

    aero_out = (
        aero.groupby(["co_muni_ori", "year", "month"], as_index=False)
        .agg(
            {
                "aero_pass": "sum",
                "aero_pass_week": "sum",
                "aero_conec": "sum",
                "co_muni_des": "nunique",
            }
        )
        .rename(
            columns={
                "co_muni_ori": "ibge_code",
                "aero_pass": "air_pass_out",
                "aero_pass_week": "air_pass_week_out",
                "aero_conec": "air_conec_out",
                "co_muni_des": "air_destinations_n",
            }
        )
    )

    aero_in = (
        aero.groupby(["co_muni_des", "year", "month"], as_index=False)
        .agg(
            {
                "aero_pass": "sum",
                "aero_pass_week": "sum",
                "aero_conec": "sum",
                "co_muni_ori": "nunique",
            }
        )
        .rename(
            columns={
                "co_muni_des": "ibge_code",
                "aero_pass": "air_pass_in",
                "aero_pass_week": "air_pass_week_in",
                "aero_conec": "air_conec_in",
                "co_muni_ori": "air_origins_n",
            }
        )
    )

    return aero_out.merge(
        aero_in,
        on=["ibge_code", "year", "month"],
        how="outer",
    ).fillna(0)


def build_fluvi_road_features(valid_ibge) -> pd.DataFrame:
    fluvi = clean_columns(pd.read_parquet(FLUVI_FILE))
    numeric_cols = [
        "co_muni_ori",
        "co_muni_des",
        "fluv_conec",
        "road_conec",
        "tot_conec",
        "irregular_conec",
    ]
    for col in numeric_cols:
        fluvi[col] = pd.to_numeric(fluvi[col], errors="coerce")

    fluvi = fluvi.dropna(subset=["co_muni_ori", "co_muni_des"]).copy()
    fluvi["co_muni_ori"] = fluvi["co_muni_ori"].astype(int)
    fluvi["co_muni_des"] = fluvi["co_muni_des"].astype(int)
    fluvi = fluvi[
        fluvi["co_muni_ori"].isin(valid_ibge) | fluvi["co_muni_des"].isin(valid_ibge)
    ].copy()

    fluvi_out = (
        fluvi.groupby("co_muni_ori", as_index=False)
        .agg(
            {
                "fluv_conec": "sum",
                "road_conec": "sum",
                "tot_conec": "sum",
                "irregular_conec": "sum",
                "co_muni_des": "nunique",
            }
        )
        .rename(
            columns={
                "co_muni_ori": "ibge_code",
                "fluv_conec": "fluv_conec_out",
                "road_conec": "road_conec_out",
                "tot_conec": "tot_conec_out",
                "irregular_conec": "irregular_conec_out",
                "co_muni_des": "network_destinations_n",
            }
        )
    )

    fluvi_in = (
        fluvi.groupby("co_muni_des", as_index=False)
        .agg(
            {
                "fluv_conec": "sum",
                "road_conec": "sum",
                "tot_conec": "sum",
                "irregular_conec": "sum",
                "co_muni_ori": "nunique",
            }
        )
        .rename(
            columns={
                "co_muni_des": "ibge_code",
                "fluv_conec": "fluv_conec_in",
                "road_conec": "road_conec_in",
                "tot_conec": "tot_conec_in",
                "irregular_conec": "irregular_conec_in",
                "co_muni_ori": "network_origins_n",
            }
        )
    )

    return fluvi_out.merge(fluvi_in, on="ibge_code", how="outer").fillna(0)


def build_model_dataframe() -> pd.DataFrame:
    df = clean_columns(pd.read_csv(COMBINED_FILE))
    df["municipio"] = df["municipio"].apply(normalize_name)

    for col in ["year", "week", "cases"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    covariate_cols = ["humidity", "idhm", "population", "rainfall", "temperature"]
    for col in covariate_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["municipio", "year", "week", "cases"]).copy()
    df["year"] = df["year"].astype(int)
    df["week"] = df["week"].astype(int)
    df = df[(df["year"] >= DATA_START_YEAR) & (df["year"] <= DATA_END_YEAR)].copy()

    lookup = build_municipio_lookup()
    hub = build_hub_features()

    df = df.merge(lookup, on="municipio", how="left", validate="many_to_one")
    missing_lookup = sorted(df[df["ibge_code"].isna()]["municipio"].unique().tolist())
    if missing_lookup:
        print("Municipios missing IBGE lookup:", missing_lookup[:20])

    df = df.dropna(subset=["ibge_code"]).copy()
    df["ibge_code"] = df["ibge_code"].astype(int)
    valid_ibge = set(df["ibge_code"].unique())

    df["date"] = iso_week_to_date(df["year"], df["week"])
    df = df.dropna(subset=["date"]).copy()
    df["month"] = df["date"].dt.month.astype(int)

    df = df.merge(
        hub.drop(columns=["municipio"]),
        on="ibge_code",
        how="left",
        validate="many_to_one",
    )

    aero = build_aero_features(valid_ibge)
    fluvi = build_fluvi_road_features(valid_ibge)

    df = df.merge(
        aero,
        on=["ibge_code", "year", "month"],
        how="left",
        validate="many_to_one",
    ).merge(
        fluvi,
        on="ibge_code",
        how="left",
        validate="many_to_one",
    )

    fill_zero_cols = [
        "air_pass_in",
        "air_pass_out",
        "air_pass_week_in",
        "air_pass_week_out",
        "air_conec_in",
        "air_conec_out",
        "air_origins_n",
        "air_destinations_n",
        "road_conec_in",
        "road_conec_out",
        "fluv_conec_in",
        "fluv_conec_out",
        "tot_conec_in",
        "tot_conec_out",
        "irregular_conec_in",
        "irregular_conec_out",
        "network_origins_n",
        "network_destinations_n",
    ]
    for col in fill_zero_cols:
        if col in df.columns:
            df[col] = df[col].fillna(0)

    df["population"] = df["population"].fillna(df["population_2022"])
    df["density_2022"] = df["density_2022"].fillna(0)
    df["hub_degree"] = df["hub_degree"].fillna(0)
    df["hub_proximity_pct"] = df["hub_proximity_pct"].fillna(0)

    df = df.sort_values(["municipio", "date"]).reset_index(drop=True)
    df["cases_lag"] = df.groupby("municipio")["cases"].shift(CASE_LAG_WEEKS)
    df["log_cases_lag"] = np.log1p(df["cases_lag"])
    df["log_population"] = np.log1p(df["population"])
    df["log_density_2022"] = np.log1p(df["density_2022"])

    for col in SKEWED_FEATURES:
        df[col] = np.log1p(df[col].clip(lower=0))

    required = ["cases", *BASE_COVARIATES]
    df = df.dropna(subset=required).copy()
    df["cases"] = df["cases"].clip(lower=0).astype(int)

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
        pm.NegativeBinomial("cases", mu=mu, alpha=alpha_nb, observed=inputs["y"], dims="obs_id")

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
    model, trace, posterior_predictive = fit_model(inputs)

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
    print(az.summary(trace, var_names=["intercept", "beta", "sigma_municipio", "sigma_week", "sigma_year", "alpha_nb"]))

    if SAVE_OUTPUTS:
        metrics_path = os.path.join(BASE_DIR, "base_model_metrics.csv")
        predictions_path = os.path.join(BASE_DIR, "base_model_predictions.csv")

        pd.DataFrame([{**metrics, **diagnostics}]).to_csv(metrics_path, index=False)

        output_df = inputs["df"][
            ["municipio", "ibge_code", "year", "week", "date", "cases"]
        ].copy()
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
