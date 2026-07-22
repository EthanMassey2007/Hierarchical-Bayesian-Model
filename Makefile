PYTHON ?= python3
RSCRIPT ?= Rscript

.PHONY: check python-check r-check setup-r m0 m5 s1 s2 s6 s7 diagnostics maps

check: python-check r-check

python-check:
	PYTHONPYCACHEPREFIX=/tmp/hbm_pycache_check $(PYTHON) -m py_compile base_model/*.py unexplained_effects_map.py

r-check:
	$(RSCRIPT) -e 'files <- list.files("spatial_R", pattern = "[.]R$$", full.names = TRUE); invisible(lapply(files, parse))'

setup-r:
	$(RSCRIPT) spatial_R/install_packages.R

m0:
	$(PYTHON) base_model/base_model.py

m5:
	$(PYTHON) base_model/base_model_lag_cases_weather.py

s1:
	$(RSCRIPT) spatial_R/spatial_inla_model_s1.R

s2:
	$(RSCRIPT) spatial_R/spatial_inla_model_s2.R

s6:
	$(RSCRIPT) spatial_R/spatial_inla_model_s6_rainfall_region.R

s7:
	$(RSCRIPT) spatial_R/spatial_inla_model_s7_temperature_region.R

diagnostics:
	$(RSCRIPT) spatial_R/s2_morans_i_diagnostic.R

maps:
	$(RSCRIPT) spatial_R/map_s2_unexplained_effects.R
	$(RSCRIPT) spatial_R/map_s6_rainfall_region_effect.R
	$(RSCRIPT) spatial_R/map_s7_temperature_region_effect.R
