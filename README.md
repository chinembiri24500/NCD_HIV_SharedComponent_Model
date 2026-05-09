# SADC HIV-NCD Bayesian spatial analysis

This repository contains the R code used for the revised manuscript:

**Uncovering Spatial Co-dependence and Potential Syndemic Overlap: A Bayesian Stacked-Likelihood Model of HIV Prevalence and NCD Mortality in Southern Africa**

The analysis estimates a primary stacked-likelihood Bayesian areal model and a shared-component sensitivity model for HIV prevalence and major NCD mortality in 13 Southern African Development Community countries, 2000-2019.

## Repository structure

```text
.
├── sadc_hiv_ncd_inla_analysis.R
├── README.md
├── data/
│   ├── sadc2025.csv
│   ├── HIV_df.csv
│   ├── sadc.shp
│   ├── sadc.dbf
│   ├── sadc.shx
│   └── sadc.prj
└── outputs/                 # created by the script
    ├── figures/
    ├── tables/
    └── models/
```

## Required input files

The script expects the following files in a `data/` folder:

1. `sadc2025.csv`: NCD panel data. Required columns are:
   - `country`
   - `year`
   - `observed_deaths`
   - `expected_deaths`
   - `ncdstype`
   - `gdpercpta`
   - `urbrte`
   - `hcexp`
   - `avrgprecp`
   - `avrgtemp`

2. `HIV_df.csv`: HIV prevalence data. Required columns are:
   - `country`
   - `year`
   - `hiv_prev`

3. `sadc.shp` and associated shapefile sidecar files. The shapefile must contain a `country` column.

## R version and packages

The script requires R >= 4.1 because it uses the base R pipe operator (`|>`).

Required R packages:

- `INLA`
- `sf`
- `spdep`
- `dplyr`
- `ggplot2`
- `readr`
- `stringr`
- `tmap`
- `scales`
- `tidyr`
- `purrr`
- `tibble`
- `data.table`
- `grid`

INLA is not installed from CRAN. It can be installed using:

```r
install.packages(
  "INLA",
  repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
  dep = TRUE
)
```

## How to run

Place the input files in the `data/` folder, set the working directory to the repository root, and run:

```r
source("sadc_hiv_ncd_inla_analysis.R")
```

The script creates the `outputs/` folder automatically.

## Main analytical steps

The script follows the revised manuscript structure:

1. Data preparation and harmonisation
2. Distance-based spatial adjacency construction
3. Primary stacked-likelihood Poisson-Beta model
4. Model-fit table for primary M1-M4 tiers
5. HIV fixed-effect table
6. Shared-component sensitivity model for M1 and M4
7. Shared-component fit and loading-parameter tables
8. Shared spatial component map
9. Shared temporal component figure
10. Residual co-clustering map after shared-component adjustment
11. HIV and NCD observed-versus-fitted calibration diagnostics

## Main outputs

### Tables

- `outputs/tables/table_1_primary_model_fit.csv`
- `outputs/tables/table_2_primary_hyperparameters_raw.csv`
- `outputs/tables/table_3_shared_component_sensitivity.csv`
- `outputs/tables/table_4_hiv_fixed_effects.csv`

### Figures

- `outputs/figures/figure_2_shared_spatial_component.png`
- `outputs/figures/figure_3_shared_temporal_component.png`
- `outputs/figures/figure_4_residual_coclustering_shared_component.png`
- `outputs/figures/figure_5_hiv_observed_fitted_shared_component.png`
- `outputs/figures/figure_6_ncd_observed_fitted_shared_component.png`

### Model objects

- `outputs/models/primary_stacked_models.rds`
- `outputs/models/primary_M1.rds`
- `outputs/models/primary_M4.rds`
- `outputs/models/shared_component_M1.rds`
- `outputs/models/shared_component_M4.rds`

## Notes on interpretation

The primary stacked-likelihood model estimates outcome-specific spatial and temporal latent fields for HIV prevalence and NCD mortality. The shared-component sensitivity model evaluates whether the two outcomes load onto common BYM2 spatial and RW1 temporal fields. The outputs should be interpreted as evidence of spatial co-patterning and potential syndemic overlap, not as proof of individual-level comorbidity or causal biological interaction.

## Suggested citation

After archiving the repository on Zenodo, replace this section with the final repository DOI.
