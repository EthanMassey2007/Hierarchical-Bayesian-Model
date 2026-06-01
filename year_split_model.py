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
# Helpers
# =========================================================
def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = df.columns.str.strip().str.lower()
    return df


def normalize_municipio_name(name):
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


def weighted_average(values, weights):
    return np.average(values, weights=weights)


def da_mean(x):
    return float(np.asarray(x).mean())


def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)

    mae = float(np.mean(np.abs(y_true - y_pred)))
    rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))

    total_actual = float(np.sum(np.abs(y_true)))
    wape = float(np.sum(np.abs(y_true - y_pred)) / max(total_actual, 1e-9))
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


def summarize_diagnostics(trace):
    rhat = az.rhat(trace)
    ess_bulk = az.ess(trace, method="bulk")
    ess_tail = az.ess(trace, method="tail")

    rhat_values = [
        da_mean(rhat["beta"]),
        da_mean(rhat["alpha_group"]),
        da_mean(rhat["sigma_group"]),
    ]
    ess_bulk_values = [
        da_mean(ess_bulk["beta"]),
        da_mean(ess_bulk["alpha_group"]),
        da_mean(ess_bulk["sigma_group"]),
    ]
    ess_tail_values = [
        da_mean(ess_tail["beta"]),
        da_mean(ess_tail["alpha_group"]),
        da_mean(ess_tail["sigma_group"]),
    ]
    weights = [2, 2, 1]

    return {
        "weighted_rhat": weighted_average(rhat_values, weights),
        "weighted_ess_bulk": weighted_average(ess_bulk_values, weights),
        "weighted_ess_tail": weighted_average(ess_tail_values, weights),
    }


# =========================================================
# Config
# =========================================================
TARGET_STATE = "RJ"

TRAIN_START_YEAR = 2017
TRAIN_END_YEAR = 2022
TEST_START_YEAR = 2023
TEST_END_YEAR = 2025

DATA_START_YEAR = TRAIN_START_YEAR
DATA_END_YEAR = TEST_END_YEAR

START_LAG = 1
END_LAG = 1
LAG_VALUES = list(range(START_LAG, END_LAG + 1))

MAKE_PLOTS = True
APPLY_LOG1P_TO_SKEWED_FEATURES = True
USE_FULL_COVARIATE_SET = True
SAVE_RESULTS_CSV = True

USE_AIR_FEATURES = False

ACCURACY_LABEL = "WAPE-based accuracy (%)"

DRAWS = 1000
TUNE = 2000
CHAINS = 4
CORES = 4
TARGET_ACCEPT = 0.98
RANDOM_SEED = 42

SKEWED_FEATURES = [
    "road_conec_in",
    "fluv_conec_in",
]

BASE_COVARIATES_NO_LAG = [
    "rainfall",
    "humidity",
    "temperature",
    "idhm",
    "time_idx",
]

FULL_COVARIATES_NO_LAG = [
    "rainfall",
    "humidity",
    "temperature",
    "idhm",
    "road_conec_in",
    "fluv_conec_in",
    "time_idx",
]


# =========================================================
# File paths
# =========================================================
base_dir = os.path.dirname(__file__)
data_dir = os.path.join(base_dir, "data")

combined_file = os.path.join(data_dir, "complete_combined_datasets.csv")

municipios_file = os.path.join(data_dir, "municipios.csv")
aero_file = os.path.join(data_dir, "aero_anac_2017_2023.parquet")
fluvi_file = os.path.join(data_dir, "fluvi_road_ibge.parquet")
hub_file = os.path.join(data_dir, "hub_pop_density.csv")


# =========================================================
# Shared preprocessing
# =========================================================
def restrict_years(d: pd.DataFrame) -> pd.DataFrame:
    d = d.dropna(subset=["year", "week"]).copy()
    d = d[(d["year"] >= DATA_START_YEAR) & (d["year"] <= DATA_END_YEAR)].copy()
    return d


def build_base_dataframe() -> pd.DataFrame:
    combined_df = clean_columns(pd.read_csv(combined_file))
    municipios_df = clean_columns(pd.read_csv(municipios_file))
    hub_df = clean_columns(pd.read_csv(hub_file))
    fluvi_df = clean_columns(pd.read_parquet(fluvi_file))

    aero_df = None
    if USE_AIR_FEATURES:
        aero_df = clean_columns(pd.read_parquet(aero_file))

    combined_df["municipio"] = combined_df["municipio"].apply(normalize_municipio_name)

    municipios_df["state_code"] = municipios_df["city"].astype(str).str.extract(
        r"/([A-Z]{2})$",
        expand=False,
    )
    municipios_df["city_clean"] = municipios_df["city"].apply(normalize_municipio_name)
    municipios_df["ibgeid"] = pd.to_numeric(
        municipios_df["ibgeid"], errors="coerce"
    ).astype("Int64")
    municipios_df = municipios_df[municipios_df["state_code"] == TARGET_STATE].copy()

    lookup_df = (
        municipios_df[["city_clean", "ibgeid"]]
        .dropna()
        .drop_duplicates()
        .rename(columns={"city_clean": "municipio", "ibgeid": "ibge_code"})
    )

    hub_df = hub_df[hub_df["uf"].astype(str).str.upper() == TARGET_STATE].copy()
    hub_df["municipio"] = hub_df["nm_municipio"].apply(normalize_municipio_name)
    hub_df["ibge_code"] = pd.to_numeric(hub_df["co_ibge"], errors="coerce")
    hub_lookup_df = (
        hub_df[["municipio", "ibge_code"]]
        .dropna()
        .drop_duplicates()
        .assign(ibge_code=lambda d: d["ibge_code"].astype(int))
    )
    lookup_df = (
        pd.concat([lookup_df, hub_lookup_df], ignore_index=True)
        .drop_duplicates(subset=["municipio"], keep="first")
        .reset_index(drop=True)
    )

    dup_counts = lookup_df.groupby("municipio")["ibge_code"].nunique()
    ambiguous_after_filter = dup_counts[dup_counts > 1]
    print("RJ municipios in lookup:", len(lookup_df))
    print("Ambiguous names after state filter:")
    print(ambiguous_after_filter)

    for col in ["year", "week", "cases", "temperature", "humidity", "rainfall", "idhm"]:
        combined_df[col] = pd.to_numeric(combined_df[col], errors="coerce")

    combined_df = restrict_years(combined_df)

    print("Core years after restriction:")
    print(sorted(combined_df["year"].dropna().astype(int).unique().tolist()))

    df = combined_df.groupby(["municipio", "year", "week"], as_index=False).agg(
        {
            "cases": "sum",
            "temperature": "mean",
            "humidity": "mean",
            "rainfall": "mean",
            "idhm": "mean",
        }
    )

    print("Merged weekly row count:", len(df))

    df = df.merge(
        lookup_df,
        on="municipio",
        how="left",
        validate="many_to_one",
    )

    missing_ibge = df[df["ibge_code"].isna()]["municipio"].drop_duplicates().sort_values()
    print("Municipios with no IBGE match:", len(missing_ibge))
    if len(missing_ibge) > 0:
        print(missing_ibge.tolist()[:50])

    df = df.dropna(
        subset=["cases", "temperature", "humidity", "rainfall", "idhm", "ibge_code"]
    ).copy()
    df["ibge_code"] = df["ibge_code"].astype(int)

    df["date"] = iso_week_to_date(df["year"], df["week"])
    df = df.dropna(subset=["date"]).copy()
    df["month"] = df["date"].dt.month.astype(int)

    print("Weekly years used:", sorted(df["year"].dropna().astype(int).unique().tolist()))

    valid_ibge = set(df["ibge_code"].unique())

    if USE_AIR_FEATURES:
        aero_df = aero_df.rename(columns={"ano": "year", "mes": "month"})
        for col in [
            "year",
            "month",
            "co_muni_ori",
            "co_muni_des",
            "aero_pass",
            "aero_pass_week",
            "aero_conec",
        ]:
            aero_df[col] = pd.to_numeric(aero_df[col], errors="coerce")

        aero_df = aero_df.dropna(
            subset=["year", "month", "co_muni_ori", "co_muni_des"]
        ).copy()
        aero_df["year"] = aero_df["year"].astype(int)
        aero_df["month"] = aero_df["month"].astype(int)
        aero_df["co_muni_ori"] = aero_df["co_muni_ori"].astype(int)
        aero_df["co_muni_des"] = aero_df["co_muni_des"].astype(int)
        aero_df = aero_df[
            (aero_df["year"] >= DATA_START_YEAR) & (aero_df["year"] <= TRAIN_END_YEAR)
        ].copy()

        print("Air years used:", sorted(aero_df["year"].unique().tolist()))

        aero_out = (
            aero_df.groupby(["co_muni_ori", "year", "month"], as_index=False)
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
            aero_df.groupby(["co_muni_des", "year", "month"], as_index=False)
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

        aero_features = aero_out.merge(
            aero_in,
            on=["ibge_code", "year", "month"],
            how="outer",
        ).fillna(0)
        aero_features = aero_features[aero_features["ibge_code"].isin(valid_ibge)].copy()

        df = df.merge(
            aero_features,
            on=["ibge_code", "year", "month"],
            how="left",
            validate="many_to_one",
        )

    for col in [
        "co_muni_ori",
        "co_muni_des",
        "fluv_conec",
        "road_conec",
        "tot_conec",
        "irregular_conec",
    ]:
        fluvi_df[col] = pd.to_numeric(fluvi_df[col], errors="coerce")

    fluvi_df = fluvi_df.dropna(subset=["co_muni_ori", "co_muni_des"]).copy()
    fluvi_df["co_muni_ori"] = fluvi_df["co_muni_ori"].astype(int)
    fluvi_df["co_muni_des"] = fluvi_df["co_muni_des"].astype(int)

    fluvi_out = (
        fluvi_df.groupby("co_muni_ori", as_index=False)
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
        fluvi_df.groupby("co_muni_des", as_index=False)
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

    fluvi_features = fluvi_out.merge(
        fluvi_in,
        on="ibge_code",
        how="outer",
    ).fillna(0)
    fluvi_features = fluvi_features[fluvi_features["ibge_code"].isin(valid_ibge)].copy()

    df = df.merge(
        fluvi_features,
        on="ibge_code",
        how="left",
        validate="many_to_one",
    )

    spatial_cols = [
        "air_pass_out",
        "air_pass_week_out",
        "air_conec_out",
        "air_destinations_n",
        "air_pass_in",
        "air_pass_week_in",
        "air_conec_in",
        "air_origins_n",
        "fluv_conec_out",
        "road_conec_out",
        "tot_conec_out",
        "irregular_conec_out",
        "network_destinations_n",
        "fluv_conec_in",
        "road_conec_in",
        "tot_conec_in",
        "irregular_conec_in",
        "network_origins_n",
    ]
    for col in spatial_cols:
        if col in df.columns:
            df[col] = df[col].fillna(0)

    print("Row count after spatial merges:", len(df))

    df = df.sort_values(["municipio", "date"]).reset_index(drop=True)
    return df


# =========================================================
# Lag-specific prep with train/test split
# =========================================================
def prepare_lagged_split_dataframe(df_base: pd.DataFrame, lag_weeks: int):
    df = df_base.copy()

    lag_cases_col = f"cases_lag{lag_weeks}"
    lag_log_col = f"log_cases_lag{lag_weeks}"

    df[lag_cases_col] = df.groupby("municipio")["cases"].shift(lag_weeks)
    df[lag_log_col] = np.log1p(df[lag_cases_col])
    df = df.dropna(subset=[lag_cases_col]).copy()

    df["time_idx"] = ((df["date"] - df["date"].min()).dt.days.astype(float) / 7.0)
    df["week_of_year"] = df["week"].astype(int)

    week_levels = sorted(df["week_of_year"].unique().tolist())
    week_to_idx = {w: i for i, w in enumerate(week_levels)}
    df["week_idx"] = df["week_of_year"].map(week_to_idx).astype(int)

    covariates_no_lag = (
        FULL_COVARIATES_NO_LAG if USE_FULL_COVARIATE_SET else BASE_COVARIATES_NO_LAG
    )
    covariates = covariates_no_lag + [lag_log_col]

    missing_covs = [c for c in covariates if c not in df.columns]
    if missing_covs:
        raise ValueError(f"Missing covariates: {missing_covs}")

    if APPLY_LOG1P_TO_SKEWED_FEATURES:
        for col in SKEWED_FEATURES:
            if col in df.columns:
                df[col] = np.log1p(df[col])

    train_df = df[
        (df["year"] >= TRAIN_START_YEAR) & (df["year"] <= TRAIN_END_YEAR)
    ].copy()
    test_df = df[
        (df["year"] >= TEST_START_YEAR) & (df["year"] <= TEST_END_YEAR)
    ].copy()

    if train_df.empty:
        raise ValueError("Training split is empty after preprocessing.")
    if test_df.empty:
        raise ValueError("Testing split is empty after preprocessing.")

    train_municipios = np.sort(train_df["municipio"].unique())
    municipio_to_idx = {m: i for i, m in enumerate(train_municipios)}

    train_df["municipio_idx"] = train_df["municipio"].map(municipio_to_idx).astype(int)
    test_df = test_df[test_df["municipio"].isin(municipio_to_idx)].copy()

    if test_df.empty:
        raise ValueError(
            "Testing split became empty after restricting to training municipios."
        )

    test_df["municipio_idx"] = test_df["municipio"].map(municipio_to_idx).astype(int)

    scaler = StandardScaler()
    train_df.loc[:, covariates] = scaler.fit_transform(train_df[covariates])
    test_df.loc[:, covariates] = scaler.transform(test_df[covariates])

    print(f"Using case lag of {lag_weeks} week(s)")
    print("Train rows:", len(train_df))
    print("Test rows:", len(test_df))
    print("Train years:", sorted(train_df["year"].astype(int).unique().tolist()))
    print("Test years:", sorted(test_df["year"].astype(int).unique().tolist()))

    print("Train covariate correlation matrix:")
    print(train_df[covariates].corr())

    for label, frame in [("train", train_df), ("test", test_df)]:
        y_values = frame["cases"].to_numpy(dtype=np.int64)
        zero_pct = np.mean(y_values == 0) * 100
        print(f"{label.title()} zero percentage: {zero_pct:.2f}%")
        print(
            f"{label.title()} mean/variance/min/max:",
            float(frame["cases"].mean()),
            float(frame["cases"].var()),
            int(frame["cases"].min()),
            int(frame["cases"].max()),
        )

    return {
        "train_df": train_df,
        "test_df": test_df,
        "covariates": covariates,
        "lag_log_col": lag_log_col,
        "municipios": train_municipios,
        "week_levels": week_levels,
    }


def build_model_inputs(df: pd.DataFrame, covariates: list[str]) -> dict:
    return {
        "y": df["cases"].to_numpy(dtype=np.int64),
        "X": df[covariates].to_numpy(dtype=np.float64),
        "group_idx": df["municipio_idx"].to_numpy(dtype=np.int32),
        "week_idx": df["week_idx"].to_numpy(dtype=np.int32),
    }


# =========================================================
# Model fit for one lag
# =========================================================
def fit_one_lag(df_base: pd.DataFrame, lag_weeks: int):
    prepared = prepare_lagged_split_dataframe(df_base, lag_weeks)

    train_df = prepared["train_df"]
    test_df = prepared["test_df"]
    covariates = prepared["covariates"]
    lag_log_col = prepared["lag_log_col"]
    municipios = prepared["municipios"]
    week_levels = prepared["week_levels"]

    train_inputs = build_model_inputs(train_df, covariates)
    test_inputs = build_model_inputs(test_df, covariates)

    coords = {
        "municipio": municipios,
        "covariate": covariates,
        "week_level": np.array(week_levels),
        "obs_id": np.arange(len(train_df)),
    }

    lag_col_idx = covariates.index(lag_log_col)

    print("Train X shape:", train_inputs["X"].shape)
    print("Train y shape:", train_inputs["y"].shape)
    print("Test X shape:", test_inputs["X"].shape)
    print("Test y shape:", test_inputs["y"].shape)
    print("Lag covariate:", lag_log_col, "at column index", lag_col_idx)

    with pm.Model(coords=coords) as model:
        X_data = pm.Data("X_data", train_inputs["X"], dims=("obs_id", "covariate"))
        group_data = pm.Data("group_data", train_inputs["group_idx"], dims="obs_id")
        week_data = pm.Data("week_data", train_inputs["week_idx"], dims="obs_id")

        alpha_global = pm.Normal(
            "alpha_global",
            mu=np.log(np.maximum(train_inputs["y"].mean(), 1.0)),
            sigma=3.0,
        )
        sigma_group = pm.HalfNormal("sigma_group", sigma=3.0)

        z_group_raw = pm.Normal("z_group_raw", mu=0.0, sigma=1.0, dims="municipio")
        z_group = pm.Deterministic(
            "z_group",
            z_group_raw - pm.math.mean(z_group_raw),
            dims="municipio",
        )
        alpha_group = pm.Deterministic(
            "alpha_group",
            alpha_global + sigma_group * z_group,
            dims="municipio",
        )

        sigma_week = pm.HalfNormal("sigma_week", sigma=1.0)
        week_raw = pm.Normal("week_raw", mu=0.0, sigma=1.0, dims="week_level")
        week_effect = pm.Deterministic(
            "week_effect",
            sigma_week * (week_raw - pm.math.mean(week_raw)),
            dims="week_level",
        )

        beta = pm.Normal("beta", mu=0.0, sigma=1.5, dims="covariate")

        eta = (
            alpha_group[group_data]
            + week_effect[week_data]
            + pm.math.dot(X_data, beta)
        )
        mu = pm.Deterministic("mu", pm.math.exp(eta), dims="obs_id")

        alpha_nb = pm.Exponential("alpha_nb", lam=1.0)

        zi_intercept = pm.Normal("zi_intercept", mu=0.0, sigma=1.5)
        zi_beta_lag = pm.Normal("zi_beta_lag", mu=0.0, sigma=1.5)
        logit_psi = zi_intercept + zi_beta_lag * X_data[:, lag_col_idx]
        psi = pm.Deterministic("psi", pm.math.sigmoid(logit_psi), dims="obs_id")

        pm.ZeroInflatedNegativeBinomial(
            "y_obs",
            psi=psi,
            mu=mu,
            alpha=alpha_nb,
            observed=train_inputs["y"],
            dims="obs_id",
        )

        trace = pm.sample(
            draws=DRAWS,
            tune=TUNE,
            chains=CHAINS,
            cores=CORES,
            target_accept=TARGET_ACCEPT,
            init="jitter+adapt_diag",
            random_seed=RANDOM_SEED,
            return_inferencedata=True,
            progressbar=True,
            idata_kwargs={"log_likelihood": True},
        )

        train_ppc = pm.sample_posterior_predictive(
            trace,
            var_names=["y_obs"],
            random_seed=RANDOM_SEED,
            progressbar=True,
            return_inferencedata=False,
        )

        pm.set_data(
            {
                "X_data": test_inputs["X"],
                "group_data": test_inputs["group_idx"],
                "week_data": test_inputs["week_idx"],
            },
            coords={"obs_id": np.arange(len(test_df))},
        )

        test_ppc = pm.sample_posterior_predictive(
            trace,
            var_names=["y_obs"],
            random_seed=RANDOM_SEED,
            progressbar=True,
            return_inferencedata=False,
        )

    train_draws = np.asarray(train_ppc["y_obs"])
    test_draws = np.asarray(test_ppc["y_obs"])

    train_pred_mean = train_draws.mean(axis=(0, 1))
    train_pred_std = train_draws.std(axis=(0, 1))
    test_pred_mean = test_draws.mean(axis=(0, 1))
    test_pred_std = test_draws.std(axis=(0, 1))

    train_metrics = compute_metrics(train_inputs["y"], train_pred_mean)
    test_metrics = compute_metrics(test_inputs["y"], test_pred_mean)
    diag_summary = summarize_diagnostics(trace)

    print(f"Lag {lag_weeks} train metrics:")
    print(
        f"  accuracy_pct={train_metrics['accuracy_pct']:.2f}, "
        f"rmse={train_metrics['rmse']:.4f}, "
        f"mae={train_metrics['mae']:.4f}, "
        f"r2={train_metrics['r2']:.4f}, "
        f"wape={train_metrics['wape']:.4f}"
    )
    print(f"Lag {lag_weeks} test metrics:")
    print(
        f"  accuracy_pct={test_metrics['accuracy_pct']:.2f}, "
        f"rmse={test_metrics['rmse']:.4f}, "
        f"mae={test_metrics['mae']:.4f}, "
        f"r2={test_metrics['r2']:.4f}, "
        f"wape={test_metrics['wape']:.4f}"
    )

    summary_main = az.summary(
        trace,
        var_names=[
            "alpha_global",
            "sigma_group",
            "sigma_week",
            "alpha_nb",
            "zi_intercept",
            "zi_beta_lag",
            "beta",
        ],
        round_to=4,
    )
    print(summary_main)

    print(f"Weighted Rhat: {diag_summary['weighted_rhat']:.4f}")
    print(f"Weighted ESS Bulk: {diag_summary['weighted_ess_bulk']:.2f}")
    print(f"Weighted ESS Tail: {diag_summary['weighted_ess_tail']:.2f}")

    return {
        "lag": lag_weeks,
        "trace": trace,
        "train_df": train_df,
        "test_df": test_df,
        "train_pred_mean": train_pred_mean,
        "train_pred_std": train_pred_std,
        "test_pred_mean": test_pred_mean,
        "test_pred_std": test_pred_std,
        "train_metrics": train_metrics,
        "test_metrics": test_metrics,
        "diag_summary": diag_summary,
    }


# =========================================================
# Main loop over lags
# =========================================================
def main():
    df_base = build_base_dataframe()

    lag_results = []
    best_result = None

    for lag in LAG_VALUES:
        print("\n" + "=" * 80)
        print(f"RUNNING MODEL FOR LAG = {lag}")
        print("=" * 80)

        result = fit_one_lag(df_base, lag)
        lag_results.append(
            {
                "lag": result["lag"],
                "train_accuracy_pct": result["train_metrics"]["accuracy_pct"],
                "train_mae": result["train_metrics"]["mae"],
                "train_rmse": result["train_metrics"]["rmse"],
                "train_wape": result["train_metrics"]["wape"],
                "train_r2": result["train_metrics"]["r2"],
                "test_accuracy_pct": result["test_metrics"]["accuracy_pct"],
                "test_mae": result["test_metrics"]["mae"],
                "test_rmse": result["test_metrics"]["rmse"],
                "test_wape": result["test_metrics"]["wape"],
                "test_r2": result["test_metrics"]["r2"],
                "weighted_rhat": result["diag_summary"]["weighted_rhat"],
                "weighted_ess_bulk": result["diag_summary"]["weighted_ess_bulk"],
                "weighted_ess_tail": result["diag_summary"]["weighted_ess_tail"],
                "n_train_rows": len(result["train_df"]),
                "n_test_rows": len(result["test_df"]),
            }
        )

        if (
            best_result is None
            or result["test_metrics"]["accuracy_pct"]
            > best_result["test_metrics"]["accuracy_pct"]
        ):
            best_result = result

    results_df = pd.DataFrame(lag_results).sort_values("lag").reset_index(drop=True)

    print("\nFinal lag comparison:")
    print(results_df)

    best_idx = results_df["test_accuracy_pct"].idxmax()
    best_row = results_df.loc[best_idx]
    print(
        f"\nBest lag based on test {ACCURACY_LABEL}: "
        f"lag={int(best_row['lag'])}, accuracy={best_row['test_accuracy_pct']:.2f}%"
    )

    if SAVE_RESULTS_CSV:
        results_path = os.path.join(base_dir, "lag_sweep_train_test_results.csv")
        results_df.to_csv(results_path, index=False)
        print(f"Saved lag comparison table to: {results_path}")

    if MAKE_PLOTS and best_result is not None:
        plt.figure(figsize=(10, 6))
        plt.plot(
            results_df["lag"],
            results_df["train_accuracy_pct"],
            marker="o",
            linewidth=2,
            label="Train",
        )
        plt.plot(
            results_df["lag"],
            results_df["test_accuracy_pct"],
            marker="o",
            linewidth=2,
            label="Test",
        )
        plt.xticks(results_df["lag"])
        plt.xlabel("Lag (weeks)")
        plt.ylabel(ACCURACY_LABEL)
        plt.title("Train vs test accuracy by lag")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.show()

        plt.figure(figsize=(10, 6))
        plt.plot(results_df["lag"], results_df["test_rmse"], marker="o", linewidth=2)
        plt.xticks(results_df["lag"])
        plt.xlabel("Lag (weeks)")
        plt.ylabel("Test RMSE")
        plt.title("Test RMSE vs lag")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

        y_test = best_result["test_df"]["cases"].to_numpy(dtype=float)
        pred_test = best_result["test_pred_mean"]
        pred_test_std = best_result["test_pred_std"]

        plt.figure(figsize=(8, 8))
        plt.errorbar(
            y_test,
            pred_test,
            yerr=pred_test_std,
            fmt="o",
            alpha=0.15,
            capsize=2,
        )
        line_min = min(float(y_test.min()), float(pred_test.min()))
        line_max = max(float(y_test.max()), float(pred_test.max()))
        plt.plot([line_min, line_max], [line_min, line_max], "r--")
        plt.xlabel("Actual test cases")
        plt.ylabel("Predicted test cases")
        plt.title(f"Best test fit: lag = {best_result['lag']}")
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    main()
