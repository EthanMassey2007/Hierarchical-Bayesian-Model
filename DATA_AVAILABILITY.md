# Data Availability

The data used by this project live in `data/`. The main analysis file is a
municipality-week panel for dengue modeling in Rio de Janeiro. A small set of
support files is also included for adjacency, mobility, region lookups, and map
outputs.

## Data Files

| File | Role | Source |
| --- | --- | --- |
| `data/complete_combined_datasets.csv` | Analysis-ready weekly municipal dataset with dengue cases, humidity, temperature, rainfall, population, and IDHM. | Built from InfoDengue, HDX rainfall, IBGE population estimates, and Synapse IDHM sources listed below. |
| `data/adjacency_matrix_correct.parquet` | Municipality adjacency matrix used for BYM2 spatial structure and neighboring-case lag features. | `andrezaleite/reproducibility_transportation_hubs-early_warning_surveillance_systems`. |
| `data/aero_anac_2017_2023.parquet` | Air passenger mobility data used in the S5 sensitivity model. | `andrezaleite/reproducibility_transportation_hubs-early_warning_surveillance_systems`. |
| `data/fluvi_road_ibge.parquet` | Road and fluvial connectivity data used in the S4 sensitivity model. | `andrezaleite/reproducibility_transportation_hubs-early_warning_surveillance_systems`. |
| `data/hub_pop_density.csv` | Municipality hub, IBGE region, population, and density lookup data. | `andrezaleite/reproducibility_transportation_hubs-early_warning_surveillance_systems`. |
| `data/municipios.csv` | Municipality metadata and IBGE lookup support. | `andrezaleite/reproducibility_transportation_hubs-early_warning_surveillance_systems`. |
| `data/RJ.json` | Rio de Janeiro municipality GeoJSON boundaries used for maps. | Source should be verified before public release. |

## Primary Data Sources

- Weekly dengue cases, temperature, and humidity came from InfoDengue API calls
  at the municipal level for January 1, 2010 through October 13, 2025.
- Weekly rainfall came from the Humanitarian Data Exchange Brazil subnational
  rainfall indicators for January 1, 1981 through October 11, 2025.
- Population was taken from IBGE population estimates for 2020.
- IDHM was taken from the Synapse-hosted dengue prediction dataset for 2010.
- The adjacency, air mobility, road/fluvial mobility, hub density, and
  municipality lookup files were drawn from
  <https://github.com/andrezaleite/reproducibility_transportation_hubs-early_warning_surveillance_systems>.

## Redistribution Notes

Before making this repository public, check the redistribution terms for each
raw and derived source. If a source cannot be redistributed, keep the affected
files out of the public repository and provide the steps needed to recreate
`complete_combined_datasets.csv` from the original data providers.

The only file still missing a confirmed provenance note is `data/RJ.json`.
Before submission or public release, confirm where that GeoJSON boundary file
came from and add the source here.

## References

BibTeX entries for the data sources are provided in `references.bib`.
