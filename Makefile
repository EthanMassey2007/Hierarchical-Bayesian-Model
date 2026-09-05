PYTHON ?= python3
RSCRIPT ?= Rscript
PROJECT_DIR := $(CURDIR)
PYTHONPATH := $(PROJECT_DIR)
R_BASELINE_DIR := models/r_inla/baseline
R_SPATIAL_DIR := models/r_inla/spatial
PYMC_DIR := models/pymc
ANALYSIS_DIR := scripts/analysis

.PHONY: check python-check r-check setup-r m0 m1-r m2-r m3-r m4-r m5 m5-r s1 s2 s6 s7 lag-sensitivity all-results collect-results diagnostics maps

check: python-check r-check

python-check:
	PYTHONPYCACHEPREFIX=/tmp/hbm_pycache_check $(PYTHON) -m py_compile $(PYMC_DIR)/*.py $(ANALYSIS_DIR)/*.py scripts/figures/*.py

r-check:
	HBM_PROJECT_DIR="$(PROJECT_DIR)" $(RSCRIPT) -e 'files <- c("$(ANALYSIS_DIR)/run_all_models_collect_results.R", list.files("$(R_BASELINE_DIR)", pattern = "[.]R$$", full.names = TRUE), list.files("$(R_SPATIAL_DIR)", pattern = "[.]R$$", full.names = TRUE)); invisible(lapply(files, parse))'

setup-r:
	HBM_PROJECT_DIR="$(PROJECT_DIR)" $(RSCRIPT) $(R_SPATIAL_DIR)/install_packages.R

m0:
	HBM_PROJECT_DIR="$(PROJECT_DIR)" $(PYTHON) $(PYMC_DIR)/base_model.py

m1-r:
	HBM_PROJECT_DIR="$(PROJECT_DIR)" $(RSCRIPT) $(R_BASELINE_DIR)/base_model_r_m1_covariates.R

m2-r:
	HBM_PROJECT_DIR="$(PROJECT_DIR)" $(RSCRIPT) $(R_BASELINE_DIR)/base_model_r_m2_lag_weather.R

m3-r:
	HBM_PROJECT_DIR="$(PROJECT_DIR)" $(RSCRIPT) $(R_BASELINE_DIR)/base_model_r_m3_interpolation.R

m4-r:
	HBM_PROJECT_DIR="$(PROJECT_DIR)" $(RSCRIPT) $(R_BASELINE_DIR)/base_model_r_m4_lag_cases.R

m5:
	HBM_PROJECT_DIR="$(PROJECT_DIR)" $(PYTHON) $(PYMC_DIR)/base_model_lag_cases_weather.py

m5-r:
	HBM_PROJECT_DIR="$(PROJECT_DIR)" $(RSCRIPT) $(R_BASELINE_DIR)/base_model_r_m5_lag_weather_cases.R

s1:
	HBM_PROJECT_DIR="$(PROJECT_DIR)" $(RSCRIPT) $(R_SPATIAL_DIR)/spatial_inla_model_s1.R

s2:
	HBM_PROJECT_DIR="$(PROJECT_DIR)" $(RSCRIPT) $(R_SPATIAL_DIR)/spatial_inla_model_s2.R

s6:
	HBM_PROJECT_DIR="$(PROJECT_DIR)" $(RSCRIPT) $(R_SPATIAL_DIR)/spatial_inla_model_s6_rainfall_region.R

s7:
	HBM_PROJECT_DIR="$(PROJECT_DIR)" $(RSCRIPT) $(R_SPATIAL_DIR)/spatial_inla_model_s7_temperature_region.R

lag-sensitivity:
	HBM_PROJECT_DIR="$(PROJECT_DIR)" $(RSCRIPT) $(R_SPATIAL_DIR)/lag_sensitivity_climate_s2.R

all-results:
	HBM_PROJECT_DIR="$(PROJECT_DIR)" $(RSCRIPT) $(ANALYSIS_DIR)/run_all_models_collect_results.R

collect-results:
	HBM_PROJECT_DIR="$(PROJECT_DIR)" RUN_MODEL_SCRIPTS=0 $(RSCRIPT) $(ANALYSIS_DIR)/run_all_models_collect_results.R

diagnostics:
	HBM_PROJECT_DIR="$(PROJECT_DIR)" $(RSCRIPT) $(R_SPATIAL_DIR)/s2_morans_i_diagnostic.R

maps:
	HBM_PROJECT_DIR="$(PROJECT_DIR)" $(RSCRIPT) $(R_SPATIAL_DIR)/map_s2_unexplained_effects.R
	HBM_PROJECT_DIR="$(PROJECT_DIR)" $(RSCRIPT) $(R_SPATIAL_DIR)/map_s6_rainfall_region_effect.R
	HBM_PROJECT_DIR="$(PROJECT_DIR)" $(RSCRIPT) $(R_SPATIAL_DIR)/map_s7_temperature_region_effect.R
