PYTHON ?= python3
RSCRIPT ?= Rscript

.PHONY: check python-check r-check setup-r m0 m1-r m2-r m3-r m4-r m5 m5-r s1 s2 s6 s7 lag-sensitivity all-results collect-results diagnostics maps

check: python-check r-check

python-check:
	PYTHONPYCACHEPREFIX=/tmp/hbm_pycache_check $(PYTHON) -m py_compile base_model/*.py unexplained_effects_map.py

r-check:
	$(RSCRIPT) -e 'files <- c("run_all_models_collect_results.R", list.files(pattern = "^base_model_r_.*[.]R$$"), list.files("spatial_R", pattern = "[.]R$$", full.names = TRUE)); invisible(lapply(files, parse))'

setup-r:
	$(RSCRIPT) spatial_R/install_packages.R

m0:
	$(PYTHON) base_model/base_model.py

m1-r:
	$(RSCRIPT) base_model_r_m1_covariates.R

m2-r:
	$(RSCRIPT) base_model_r_m2_lag_weather.R

m3-r:
	$(RSCRIPT) base_model_r_m3_interpolation.R

m4-r:
	$(RSCRIPT) base_model_r_m4_lag_cases.R

m5:
	$(PYTHON) base_model/base_model_lag_cases_weather.py

m5-r:
	$(RSCRIPT) base_model_r_m5_lag_weather_cases.R

s1:
	$(RSCRIPT) spatial_R/spatial_inla_model_s1.R

s2:
	$(RSCRIPT) spatial_R/spatial_inla_model_s2.R

s6:
	$(RSCRIPT) spatial_R/spatial_inla_model_s6_rainfall_region.R

s7:
	$(RSCRIPT) spatial_R/spatial_inla_model_s7_temperature_region.R

lag-sensitivity:
	$(RSCRIPT) spatial_R/lag_sensitivity_climate_s2.R

all-results:
	$(RSCRIPT) run_all_models_collect_results.R

collect-results:
	RUN_MODEL_SCRIPTS=0 $(RSCRIPT) run_all_models_collect_results.R

diagnostics:
	$(RSCRIPT) spatial_R/s2_morans_i_diagnostic.R

maps:
	$(RSCRIPT) spatial_R/map_s2_unexplained_effects.R
	$(RSCRIPT) spatial_R/map_s6_rainfall_region_effect.R
	$(RSCRIPT) spatial_R/map_s7_temperature_region_effect.R
