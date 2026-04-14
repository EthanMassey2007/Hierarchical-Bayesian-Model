import os
import re
import unicodedata
import warnings

import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
import matplotlib.pyplot as plt

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

    # Remove trailing /UF, e.g. "São Paulo/SP" -> "São Paulo"
    name = re.sub(r"/[a-z]{2}$", "", name)

    # Remove accents
    name = unicodedata.normalize("NFKD", name)
    name = "".join(c for c in name if not unicodedata.combining(c))

    # Normalize whitespace
    name = re.sub(r"\s+", " ", name).strip()

    return name


def iso_week_to_month(year_series, week_series):
    dt = pd.to_datetime(
        year_series.astype(str)
        + "-W"
        + week_series.astype(int).astype(str).str.zfill(2)
        + "-1",
        format="%G-W%V-%u",
        errors="coerce",
    )
    return dt, dt.dt.month


def weighted_average(values, weights):
    return np.average(values, weights=weights)


def da_mean(x):
    return float(np.asarray(x).mean())


# =========================================================
# Config
# =========================================================
TARGET_STATE = "RJ"
START_YEAR = 2017
END_YEAR = 2022

MAKE_PLOTS = True
APPLY_LOG1P_TO_SKEWED_FEATURES = True

BASE_COVARIATES = [
    "rainfall",
    "humidity",
    "temperature",
    "idhm",
    "air_pass_in",
]

FULL_COVARIATES = [
    "rainfall",
    "humidity",
    "temperature",
    "idhm",
    "air_pass_in",
    "road_conec_in",
    "fluv_conec_in",
]

USE_FULL_COVARIATE_SET = True

SKEWED_FEATURES = [
    "air_pass_in",
    "road_conec_in",
    "fluv_conec_in",
]


# =========================================================
# File paths
# =========================================================
base_dir = os.path.dirname(__file__)
data_dir = os.path.join(base_dir, "data")

cases_file = os.path.join(data_dir, "cases.csv")
temperature_file = os.path.join(data_dir, "temperature.csv")
humidity_file = os.path.join(data_dir, "humidity.csv")
rainfall_file = os.path.join(data_dir, "rainfall.csv")
idhm_file = os.path.join(data_dir, "idhm.csv")

municipios_file = os.path.join(data_dir, "municipios.csv")
aero_file = os.path.join(data_dir, "aero_anac_2017_2023.parquet")
fluvi_file = os.path.join(data_dir, "fluvi_road_ibge.parquet")


# =========================================================
# Load data
# =========================================================
cases_df = clean_columns(pd.read_csv(cases_file))
temp_df = clean_columns(pd.read_csv(temperature_file))
hum_df = clean_columns(pd.read_csv(humidity_file))
rain_df = clean_columns(pd.read_csv(rainfall_file))
idhm_df = clean_columns(pd.read_csv(idhm_file))

municipios_df = clean_columns(pd.read_csv(municipios_file))
aero_df = clean_columns(pd.read_parquet(aero_file))
fluvi_df = clean_columns(pd.read_parquet(fluvi_file))


# =========================================================
# Normalize municipio names in the main weekly CSVs
# =========================================================
for d in [cases_df, temp_df, hum_df, rain_df, idhm_df]:
    d["municipio"] = d["municipio"].apply(normalize_municipio_name)


# =========================================================
# Build RJ-only city -> IBGE lookup from municipios.csv
# FIX: filter by state code, not by exact city name
# =========================================================
municipios_df["state_code"] = municipios_df["city"].astype(str).str.extract(
    r"/([A-Z]{2})$",
    expand=False,
)
municipios_df["city_clean"] = municipios_df["city"].apply(normalize_municipio_name)
municipios_df["ibgeid"] = pd.to_numeric(municipios_df["ibgeid"], errors="coerce").astype("Int64")

municipios_df = municipios_df[municipios_df["state_code"] == TARGET_STATE].copy()

lookup_df = (
    municipios_df[["city_clean", "ibgeid"]]
    .dropna()
    .drop_duplicates()
    .rename(columns={"city_clean": "municipio", "ibgeid": "ibge_code"})
)

dup_counts = lookup_df.groupby("municipio")["ibge_code"].nunique()
ambiguous_after_filter = dup_counts[dup_counts > 1]

print("RJ municipios in lookup:", len(lookup_df))
print("Ambiguous names after state filter:")
print(ambiguous_after_filter)


# =========================================================
# Numeric conversion
# =========================================================
for d in [cases_df, temp_df, hum_df, rain_df, idhm_df]:
    d["year"] = pd.to_numeric(d["year"], errors="coerce").astype("Int64")
    d["week"] = pd.to_numeric(d["week"], errors="coerce").astype("Int64")

cases_df["cases"] = pd.to_numeric(cases_df["cases"], errors="coerce")
temp_df["temperature"] = pd.to_numeric(temp_df["temperature"], errors="coerce")
hum_df["humidity"] = pd.to_numeric(hum_df["humidity"], errors="coerce")
rain_df["rainfall"] = pd.to_numeric(rain_df["rainfall"], errors="coerce")
idhm_df["idhm"] = pd.to_numeric(idhm_df["idhm"], errors="coerce")


# =========================================================
# Restrict weekly core datasets to 2017-2022
# =========================================================
for name, d in [
    ("cases", cases_df),
    ("temp", temp_df),
    ("humidity", hum_df),
    ("rain", rain_df),
    ("idhm", idhm_df),
]:
    d.dropna(subset=["year", "week"], inplace=True)
    d = d[(d["year"] >= START_YEAR) & (d["year"] <= END_YEAR)].copy()

    if name == "cases":
        cases_df = d
    elif name == "temp":
        temp_df = d
    elif name == "humidity":
        hum_df = d
    elif name == "rain":
        rain_df = d
    elif name == "idhm":
        idhm_df = d

print("Core data year coverage after restriction:")
print("cases:", sorted(cases_df["year"].dropna().astype(int).unique().tolist()))
print("temp:", sorted(temp_df["year"].dropna().astype(int).unique().tolist()))
print("humidity:", sorted(hum_df["year"].dropna().astype(int).unique().tolist()))
print("rain:", sorted(rain_df["year"].dropna().astype(int).unique().tolist()))
print("idhm:", sorted(idhm_df["year"].dropna().astype(int).unique().tolist()))


# =========================================================
# Aggregate each input BEFORE merging
# =========================================================
cases_df = (
    cases_df.groupby(["municipio", "year", "week"], as_index=False)
    .agg({"cases": "sum"})
)

temp_df = (
    temp_df.groupby(["municipio", "year", "week"], as_index=False)
    .agg({"temperature": "mean"})
)

hum_df = (
    hum_df.groupby(["municipio", "year", "week"], as_index=False)
    .agg({"humidity": "mean"})
)

rain_df = (
    rain_df.groupby(["municipio", "year", "week"], as_index=False)
    .agg({"rainfall": "mean"})
)

idhm_df = (
    idhm_df.groupby(["municipio", "year", "week"], as_index=False)
    .agg({"idhm": "mean"})
)


# =========================================================
# Merge weekly core data
# =========================================================
df = (
    cases_df
    .merge(temp_df, on=["municipio", "year", "week"], how="left", validate="one_to_one")
    .merge(hum_df, on=["municipio", "year", "week"], how="left", validate="one_to_one")
    .merge(rain_df, on=["municipio", "year", "week"], how="left", validate="one_to_one")
    .merge(idhm_df, on=["municipio", "year", "week"], how="left", validate="one_to_one")
)

print("Merged weekly df columns:", df.columns.tolist())
print("Merged weekly row count:", len(df))


# =========================================================
# Add IBGE code
# =========================================================
df = df.merge(
    lookup_df,
    on="municipio",
    how="left",
    validate="many_to_one",
)

missing_ibge = df[df["ibge_code"].isna()]["municipio"].drop_duplicates().sort_values()
print(f"Municipios with no IBGE match: {len(missing_ibge)}")
if len(missing_ibge) > 0:
    print(missing_ibge.tolist()[:50])

df = df.dropna(
    subset=["cases", "temperature", "humidity", "rainfall", "idhm", "ibge_code"]
).copy()

df["ibge_code"] = df["ibge_code"].astype(int)

print("Rows after dropping missing weekly core data / IBGE:", len(df))


# =========================================================
# Add month from ISO week
# =========================================================
df["date"], df["month"] = iso_week_to_month(df["year"], df["week"])
df = df.dropna(subset=["month"]).copy()
df["month"] = df["month"].astype(int)

print("Final weekly years present:", sorted(df["year"].dropna().astype(int).unique().tolist()))


# =========================================================
# Keep only IBGE codes that survive in core weekly data
# =========================================================
valid_ibge = set(df["ibge_code"].unique())


# =========================================================
# Prepare air travel features (time-varying, 2017-2022)
# =========================================================
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

aero_df = aero_df.dropna(subset=["year", "month", "co_muni_ori", "co_muni_des"]).copy()
aero_df["year"] = aero_df["year"].astype(int)
aero_df["month"] = aero_df["month"].astype(int)
aero_df["co_muni_ori"] = aero_df["co_muni_ori"].astype(int)
aero_df["co_muni_des"] = aero_df["co_muni_des"].astype(int)

aero_df = aero_df[
    (aero_df["year"] >= START_YEAR) & (aero_df["year"] <= END_YEAR)
].copy()

print("Air data years after restriction:", sorted(aero_df["year"].unique().tolist()))

aero_out = (
    aero_df.groupby(["co_muni_ori", "year", "month"], as_index=False)
    .agg({
        "aero_pass": "sum",
        "aero_pass_week": "sum",
        "aero_conec": "sum",
        "co_muni_des": "nunique",
    })
    .rename(columns={
        "co_muni_ori": "ibge_code",
        "aero_pass": "air_pass_out",
        "aero_pass_week": "air_pass_week_out",
        "aero_conec": "air_conec_out",
        "co_muni_des": "air_destinations_n",
    })
)

aero_in = (
    aero_df.groupby(["co_muni_des", "year", "month"], as_index=False)
    .agg({
        "aero_pass": "sum",
        "aero_pass_week": "sum",
        "aero_conec": "sum",
        "co_muni_ori": "nunique",
    })
    .rename(columns={
        "co_muni_des": "ibge_code",
        "aero_pass": "air_pass_in",
        "aero_pass_week": "air_pass_week_in",
        "aero_conec": "air_conec_in",
        "co_muni_ori": "air_origins_n",
    })
)

aero_features = aero_out.merge(
    aero_in,
    on=["ibge_code", "year", "month"],
    how="outer",
).fillna(0)

aero_features = aero_features[aero_features["ibge_code"].isin(valid_ibge)].copy()


# =========================================================
# Prepare road/fluvial features (STATIC municipio-level features)
# These are NOT merged by year. They are treated as structural connectivity.
# =========================================================
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
    .agg({
        "fluv_conec": "sum",
        "road_conec": "sum",
        "tot_conec": "sum",
        "irregular_conec": "sum",
        "co_muni_des": "nunique",
    })
    .rename(columns={
        "co_muni_ori": "ibge_code",
        "fluv_conec": "fluv_conec_out",
        "road_conec": "road_conec_out",
        "tot_conec": "tot_conec_out",
        "irregular_conec": "irregular_conec_out",
        "co_muni_des": "network_destinations_n",
    })
)

fluvi_in = (
    fluvi_df.groupby("co_muni_des", as_index=False)
    .agg({
        "fluv_conec": "sum",
        "road_conec": "sum",
        "tot_conec": "sum",
        "irregular_conec": "sum",
        "co_muni_ori": "nunique",
    })
    .rename(columns={
        "co_muni_des": "ibge_code",
        "fluv_conec": "fluv_conec_in",
        "road_conec": "road_conec_in",
        "tot_conec": "tot_conec_in",
        "irregular_conec": "irregular_conec_in",
        "co_muni_ori": "network_origins_n",
    })
)

fluvi_features = fluvi_out.merge(
    fluvi_in,
    on="ibge_code",
    how="outer",
).fillna(0)

fluvi_features = fluvi_features[fluvi_features["ibge_code"].isin(valid_ibge)].copy()


# =========================================================
# Merge spatial features
# =========================================================
df = (
    df
    .merge(aero_features, on=["ibge_code", "year", "month"], how="left", validate="many_to_one")
    .merge(fluvi_features, on="ibge_code", how="left", validate="many_to_one")
)

spatial_cols = [
    "air_pass_out", "air_pass_week_out", "air_conec_out", "air_destinations_n",
    "air_pass_in", "air_pass_week_in", "air_conec_in", "air_origins_n",
    "fluv_conec_out", "road_conec_out", "tot_conec_out", "irregular_conec_out", "network_destinations_n",
    "fluv_conec_in", "road_conec_in", "tot_conec_in", "irregular_conec_in", "network_origins_n",
]

for col in spatial_cols:
    if col in df.columns:
        df[col] = df[col].fillna(0)

print("Final modeling row count:", len(df))
print("Final year counts:")
print(df.groupby("year").size())


# =========================================================
# Encode municipios
# =========================================================
df["municipio_idx"], municipios = pd.factorize(df["municipio"], sort=True)
n_groups = df["municipio_idx"].nunique()
print(f"Number of municipios: {n_groups}")

if n_groups < 2:
    raise ValueError(
        "The model has fewer than 2 municipios after preprocessing. "
        "Check the municipio lookup / merge logic."
    )


# =========================================================
# Covariates
# =========================================================
covariates = FULL_COVARIATES if USE_FULL_COVARIATE_SET else BASE_COVARIATES
print("Using covariates:", covariates)

missing_covs = [c for c in covariates if c not in df.columns]
if missing_covs:
    raise ValueError(f"Missing covariates in dataframe: {missing_covs}")

if APPLY_LOG1P_TO_SKEWED_FEATURES:
    for col in SKEWED_FEATURES:
        if col in df.columns:
            df[col] = np.log1p(df[col])

correlation_matrix = df[covariates].corr()
print("Covariate correlation matrix:")
print(correlation_matrix)

scaler = StandardScaler()
df[covariates] = scaler.fit_transform(df[covariates])


# =========================================================
# Arrays
# =========================================================
y = df["cases"].to_numpy(dtype=np.int64)
group_idx = df["municipio_idx"].to_numpy(dtype=np.int32)
X = df[covariates].to_numpy(dtype=np.float64)

print("Mean of cases:", float(df["cases"].mean()))
print("Variance of cases:", float(df["cases"].var()))
print("Min/Max of cases:", int(df["cases"].min()), int(df["cases"].max()))


# =========================================================
# Model
# =========================================================
coords = {
    "municipio": municipios,
    "covariate": covariates,
    "obs_id": np.arange(len(df)),
}

with pm.Model(coords=coords) as model:
    X_data = pm.Data("X_data", X, dims=("obs_id", "covariate"))
    group_data = pm.Data("group_data", group_idx, dims="obs_id")

    # Hierarchical intercept
    alpha_global = pm.Normal(
        "alpha_global",
        mu=np.log(np.maximum(y.mean(), 1.0)),
        sigma=1.5,
    )

    sigma_group = pm.HalfNormal("sigma_group", sigma=1.0)

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

    # Covariate coefficients
    beta = pm.Normal("beta", mu=0.0, sigma=1.0, dims="covariate")

    # Linear predictor
    eta = alpha_group[group_data] + pm.math.dot(X_data, beta)
    mu = pm.Deterministic("mu", pm.math.exp(eta), dims="obs_id")

    # Negative binomial overdispersion
    alpha_nb = pm.Exponential("alpha_nb", lam=1.0)

    # Likelihood
    pm.NegativeBinomial(
        "y_obs",
        mu=mu,
        alpha=alpha_nb,
        observed=y,
        dims="obs_id",
    )

    trace = pm.sample(
        draws=400,
        tune=600,
        chains=4,
        cores=4,
        target_accept=0.97,
        init="jitter+adapt_diag",
        random_seed=42,
        return_inferencedata=True,
        progressbar=True,
        idata_kwargs={"log_likelihood": True},
    )

    posterior_predictive = pm.sample_posterior_predictive(
        trace,
        var_names=["mu", "y_obs"],
        random_seed=42,
        progressbar=True,
    )


# =========================================================
# Diagnostics
# =========================================================
summary_main = az.summary(
    trace,
    var_names=["alpha_global", "sigma_group", "alpha_nb", "beta"],
    round_to=4,
)
print(summary_main)

summary_all = az.summary(trace, round_to=4)
print(summary_all.head())

rhat = az.rhat(trace)
ess_bulk = az.ess(trace, method="bulk")
ess_tail = az.ess(trace, method="tail")

rhat_values = [
    da_mean(rhat["beta"]),
    da_mean(rhat["alpha_group"]),
    da_mean(rhat["sigma_group"]),
    da_mean(rhat["z_group"]),
]

ess_bulk_values = [
    da_mean(ess_bulk["beta"]),
    da_mean(ess_bulk["alpha_group"]),
    da_mean(ess_bulk["sigma_group"]),
    da_mean(ess_bulk["z_group"]),
]

ess_tail_values = [
    da_mean(ess_tail["beta"]),
    da_mean(ess_tail["alpha_group"]),
    da_mean(ess_tail["sigma_group"]),
    da_mean(ess_tail["z_group"]),
]

weights = {
    "beta": 2,
    "alpha_group": 1,
    "sigma_group": 1,
    "z_group": 3,
}

weighted_rhat = weighted_average(rhat_values, list(weights.values()))
weighted_ess_bulk = weighted_average(ess_bulk_values, list(weights.values()))
weighted_ess_tail = weighted_average(ess_tail_values, list(weights.values()))

print(f"Weighted Rhat: {weighted_rhat:.4f}")
print(f"Weighted ESS Bulk: {weighted_ess_bulk:.2f}")
print(f"Weighted ESS Tail: {weighted_ess_tail:.2f}")


# =========================================================
# WAIC and LOO
# =========================================================
waic_result = az.waic(trace)
print("WAIC Result:")
print(waic_result)
print(f"WAIC (elpd_waic): {waic_result.elpd_waic:.2f}")
print(f"p_WAIC: {waic_result.p_waic:.2f}")

loo_result = az.loo(trace)
print("\nLOO Result:")
print(loo_result)
print(f"LOO (elpd_loo): {loo_result.elpd_loo:.2f}")
print(f"p_LOO: {loo_result.p_loo:.2f}")

k_vals = loo_result.pareto_k.values
high_k_idx = np.where(k_vals > 0.7)[0]
print(f"Num high Pareto-k points (> 0.7): {len(high_k_idx)}")


# =========================================================
# Effect size summaries
# =========================================================
beta_summary = az.summary(trace, var_names=["beta"], round_to=4)
beta_means = beta_summary["mean"].to_numpy()
percent_effect = np.exp(beta_means) - 1

print("\nApprox multiplicative effect of +1 SD in each covariate:")
for name, val in zip(covariates, percent_effect):
    print(f"{name}: {val * 100:.2f}% change in expected cases")


# =========================================================
# Fitted values
# =========================================================
fitted_mu_mean = (
    posterior_predictive.posterior_predictive["mu"]
    .mean(dim=["chain", "draw"])
    .values
)

fitted_mu_std = (
    posterior_predictive.posterior_predictive["mu"]
    .std(dim=["chain", "draw"])
    .values
)

ppc_y_mean = (
    posterior_predictive.posterior_predictive["y_obs"]
    .mean(dim=["chain", "draw"])
    .values
)

ppc_y_std = (
    posterior_predictive.posterior_predictive["y_obs"]
    .std(dim=["chain", "draw"])
    .values
)

print("\nFitted mean range:", float(fitted_mu_mean.min()), "to", float(fitted_mu_mean.max()))
print("Observed range:", int(y.min()), "to", int(y.max()))
print("Std of fitted means:", float(np.std(fitted_mu_mean)))

alpha_group_mean = trace.posterior["alpha_group"].mean(dim=["chain", "draw"]).values
beta_mean = trace.posterior["beta"].mean(dim=["chain", "draw"]).values
eta_covariate = X @ beta_mean

print("Alpha_group range:", float(alpha_group_mean.min()), "to", float(alpha_group_mean.max()))
print("Covariate contribution range:", float(eta_covariate.min()), "to", float(eta_covariate.max()))
print("Mean alpha_group:", float(alpha_group_mean.mean()))
print("Mean covariate contribution:", float(eta_covariate.mean()))


# =========================================================
# Prediction helper
# If municipio is supplied and known, use that municipio's posterior intercept.
# Otherwise, fall back to alpha_global.
# =========================================================
posterior_samples = trace.posterior.stack(sample=("chain", "draw"))

beta_draws = posterior_samples["beta"].values
alpha_global_draws = posterior_samples["alpha_global"].values
alpha_group_draws = posterior_samples["alpha_group"].values

municipio_to_idx = {m: i for i, m in enumerate(municipios)}


def predict_expected_cases(new_df: pd.DataFrame) -> np.ndarray:
    new_df = new_df.copy()

    for col in covariates:
        if col not in new_df.columns:
            raise ValueError(f"Missing covariate in new_df: {col}")

    if APPLY_LOG1P_TO_SKEWED_FEATURES:
        for col in SKEWED_FEATURES:
            if col in new_df.columns:
                new_df[col] = np.log1p(new_df[col])

    X_new = scaler.transform(new_df[covariates])

    if "municipio" in new_df.columns:
        normalized_names = new_df["municipio"].apply(normalize_municipio_name)
        intercept_draws = np.zeros((len(new_df), alpha_global_draws.shape[0]))

        for i, muni_name in enumerate(normalized_names):
            if muni_name in municipio_to_idx:
                intercept_draws[i, :] = alpha_group_draws[municipio_to_idx[muni_name], :]
            else:
                intercept_draws[i, :] = alpha_global_draws
    else:
        intercept_draws = np.tile(alpha_global_draws, (len(new_df), 1))

    eta_new = intercept_draws + X_new @ beta_draws
    mu_new = np.exp(eta_new)
    return mu_new.mean(axis=1)


# =========================================================
# Optional quick sanity check
# =========================================================
example_rows = df.sample(min(5, len(df)), random_state=42).copy()
example_pred = predict_expected_cases(example_rows[["municipio"] + covariates])

print("\nExample predictions:")
for muni, actual, pred in zip(example_rows["municipio"], example_rows["cases"], example_pred):
    print(f"{muni}: actual={actual}, predicted_mean={pred:.2f}")


# =========================================================
# Plots
# =========================================================
if MAKE_PLOTS:
    posterior_means = trace.posterior["beta"].mean(dim=["chain", "draw"]).values
    posterior_stds = trace.posterior["beta"].std(dim=["chain", "draw"]).values
    num_vars = len(posterior_means)

    plt.figure(figsize=(14, 6))

    plt.subplot(1, 2, 1)
    plt.errorbar(
        range(num_vars),
        posterior_means,
        yerr=posterior_stds,
        fmt="o",
        capsize=5,
    )
    plt.xticks(range(num_vars), covariates, rotation=45, ha="right")
    plt.axhline(0, color="red", linestyle="--", label="No effect")
    plt.xlabel("Covariates")
    plt.ylabel("Posterior mean coefficient ± SD")
    plt.title("Posterior Coefficients")
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.errorbar(
        y,
        fitted_mu_mean,
        yerr=fitted_mu_std,
        fmt="o",
        alpha=0.25,
        capsize=2,
    )
    line_min = min(float(y.min()), float(fitted_mu_mean.min()))
    line_max = max(float(y.max()), float(fitted_mu_mean.max()))
    plt.plot([line_min, line_max], [line_min, line_max], "r--", label="y = fitted mean")
    plt.xlabel("Actual cases")
    plt.ylabel("Fitted expected cases")
    plt.title("Fitted Expected Cases vs Actual Cases")
    plt.legend()

    plt.tight_layout()
    plt.show()

    az.plot_trace(trace, var_names=["alpha_global", "sigma_group", "alpha_nb", "beta"])
    plt.tight_layout()
    plt.show()

    az.plot_posterior(trace, var_names=["alpha_global", "sigma_group", "alpha_nb", "beta"])
    plt.tight_layout()
    plt.show()