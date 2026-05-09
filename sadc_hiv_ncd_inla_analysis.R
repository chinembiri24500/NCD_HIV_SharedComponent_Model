# =============================================================================
# SADC HIV-NCD Bayesian spatial analysis using INLA
# Primary stacked-likelihood model and shared-component sensitivity analysis
# =============================================================================
#
# Manuscript: Uncovering Spatial Co-dependence and Potential Syndemic Overlap:
# A Bayesian Stacked-Likelihood Model of HIV Prevalence and NCD Mortality in
# Southern Africa
#
# Purpose
#   This script reproduces the main model outputs used in the revised manuscript:
#   - Primary stacked Poisson-Beta model with outcome-specific BYM2/RW1 fields
#   - Shared-component sensitivity model with shared BYM2/RW1 latent fields
#   - Model-fit tables, loading-parameter table, calibration diagnostics, and maps
#
# Expected input files
#   data/sadc2025.csv  : NCD panel data with deaths, expected deaths, covariates
#   data/HIV_df.csv    : HIV prevalence data with country, year, hiv_prev columns
#   data/sadc.shp      : SADC country boundary shapefile with a country column
#
# Main outputs
#   outputs/tables/
#   outputs/figures/
#   outputs/models/
#
# Software note
#   The script uses base R pipes (|>) and therefore requires R >= 4.1.
#   INLA is not hosted on CRAN and should be installed from the INLA repository.
# =============================================================================

# =============================================================================
# 0. Configuration and packages
# =============================================================================

panel <- read_csv("data/sadc2025.csv")
hiv <- read_csv("data/HIV_df.csv")
shp <- st_read("data/sadc.shp")

set.seed(42)

base_dir <- getwd()
data_dir <- file.path(base_dir, "data")
out_dir <- file.path(base_dir, "outputs")
table_dir <- file.path(out_dir, "tables")
figure_dir <- file.path(out_dir, "figures")
model_dir <- file.path(out_dir, "models")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

required_packages <- c(
  "INLA", "sf", "spdep", "dplyr", "ggplot2", "readr", "stringr",
  "tmap", "scales", "tidyr", "purrr", "tibble", "data.table", "grid"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "The following packages are required but not installed: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(INLA)
  library(sf)
  library(spdep)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tmap)
  library(scales)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(data.table)
  library(grid)
})

# =============================================================================
# 1. Data preparation
# =============================================================================

normalize_drc <- function(x) {
  dplyr::case_when(
    x %in% c(
      "Democratic Republic of the Congo", "Dem. Rep. Congo",
      "Congo (Kinshasa)", "Congo, Dem. Rep.", "Congo Democratic Republic",
      "DR Congo"
    ) ~ "DRC",
    TRUE ~ as.character(x)
  )
}

keep_countries <- c(
  "Angola", "Botswana", "DRC", "Eswatini", "Lesotho", "Madagascar",
  "Malawi", "Mozambique", "Namibia", "South Africa", "Tanzania",
  "Zambia", "Zimbabwe"
)

panel_path <- file.path(data_dir, "sadc2025.csv")
hiv_path <- file.path(data_dir, "HIV_df.csv")
shape_path <- file.path(data_dir, "sadc.shp")

if (!file.exists(panel_path)) stop("Missing file: ", panel_path, call. = FALSE)
if (!file.exists(hiv_path)) stop("Missing file: ", hiv_path, call. = FALSE)
if (!file.exists(shape_path)) stop("Missing file: ", shape_path, call. = FALSE)

panel <- readr::read_csv(panel_path, show_col_types = FALSE)
hiv <- readr::read_csv(hiv_path, show_col_types = FALSE)

names(panel) <- tolower(names(panel))
names(hiv) <- tolower(names(hiv))

required_panel_cols <- c(
  "country", "year", "observed_deaths", "expected_deaths", "ncdstype",
  "gdpercpta", "urbrte", "hcexp", "avrgprecp", "avrgtemp"
)
required_hiv_cols <- c("country", "year", "hiv_prev")

stopifnot(all(required_panel_cols %in% names(panel)))
stopifnot(all(required_hiv_cols %in% names(hiv)))

panel <- panel |>
  dplyr::mutate(country = normalize_drc(country)) |>
  dplyr::filter(year >= 2000, year <= 2019)

hiv <- hiv |>
  dplyr::mutate(country = normalize_drc(country)) |>
  dplyr::filter(year >= 2000, year <= 2019)

dat <- panel |>
  dplyr::left_join(hiv |> dplyr::select(country, year, hiv_prev), by = c("country", "year")) |>
  tidyr::drop_na(hiv_prev) |>
  dplyr::filter(country %in% keep_countries)

scale_vars <- c("gdpercpta", "avrgtemp", "avrgprecp", "urbrte", "hcexp")
for (v in scale_vars) {
  dat[[paste0(v, "_z")]] <- as.numeric(scale(dat[[v]]))
}

dat <- dat |>
  dplyr::mutate(
    log_expected = log(expected_deaths),
    hiv_prev_p = pmin(pmax(hiv_prev / 100, 1e-4), 1 - 1e-4),
    ncdstype = factor(ncdstype)
  )

# =============================================================================
# 2. Spatial data and adjacency graph
# =============================================================================

shp <- sf::st_read(shape_path, quiet = TRUE) |>
  sf::st_make_valid()

if (!"country" %in% names(shp)) {
  stop("The shapefile must include a 'country' column.", call. = FALSE)
}

shp <- shp |>
  dplyr::mutate(country = normalize_drc(country)) |>
  dplyr::filter(country %in% keep_countries) |>
  dplyr::arrange(country)

# Distance-based neighbourhoods are used to retain non-contiguous geographies.
# This keeps Madagascar in the spatial graph as a connected node based on
# nearest regional proximity rather than as an isolated area.
centroids <- sf::st_centroid(shp)
coords <- sf::st_coordinates(centroids)
neighbours <- spdep::dnearneigh(coords, d1 = 0, d2 = 700000)

graph_path <- file.path(out_dir, "adjacency_dist.graph")
INLA::nb2INLA(graph_path, nb = neighbours)
graph <- INLA::inla.read.graph(graph_path)

dat <- dat |>
  dplyr::arrange(country, year)

shp <- shp |>
  dplyr::arrange(country) |>
  dplyr::mutate(region_id = as.integer(factor(country, levels = country)))

dat <- dat |>
  dplyr::mutate(
    region_id = as.integer(factor(country, levels = shp$country)),
    year_id = year - 2000 + 1,
    space_time_id = as.integer(factor(interaction(region_id, year_id, drop = TRUE)))
  )

stopifnot(!any(is.na(dat$region_id)))

# =============================================================================
# 3. Model priors and covariate tiers
# =============================================================================

bym2_hyper <- list(
  prec = list(prior = "pc.prec", param = c(1.5, 0.01)),
  phi = list(prior = "pc", param = c(0.5, 0.5))
)

rw1_prec <- list(
  prec = list(prior = "pc.prec", param = c(0.5, 0.01))
)

iid_prec <- list(
  prec = list(prior = "pc.prec", param = c(1, 0.01))
)

M1_vars <- c("gdpercpta_z")
M2_vars <- c("gdpercpta_z", "avrgtemp_z", "avrgprecp_z")
M3_vars <- c(M2_vars, "urbrte_z")
M4_vars <- c(M3_vars, "hcexp_z")
model_tiers <- list(M1 = M1_vars, M2 = M2_vars, M3 = M3_vars, M4 = M4_vars)

# =============================================================================
# 4. Primary stacked-likelihood model
# =============================================================================

make_interaction_columns <- function(df, covariates, factor_column) {
  if (length(covariates) == 0) return(data.frame())

  levels_factor <- levels(df[[factor_column]])
  output <- list()

  for (covariate in covariates) {
    for (level in levels_factor) {
      column_name <- paste(covariate, level, sep = "__")
      output[[column_name]] <- ifelse(df[[factor_column]] == level, df[[covariate]], 0)
    }
  }

  as.data.frame(output)
}

build_primary_stacked_model <- function(dat, graph, ncd_covariates, hiv_covariates) {
  n <- nrow(dat)

  ncd_interactions <- make_interaction_columns(dat, ncd_covariates, "ncdstype")

  effects_ncd <- c(
    list(
      intercept_ncd = rep(1, n),
      ncdstype = dat$ncdstype,
      region_id_ncd = dat$region_id,
      year_id_ncd = dat$year_id,
      space_time_ncd = dat$space_time_id
    ),
    as.list(ncd_interactions)
  )

  stack_ncd <- INLA::inla.stack(
    data = list(
      response = cbind(dat$observed_deaths, NA_real_),
      log_expected = dat$log_expected,
      country = dat$country,
      year = dat$year,
      ncdstype_diag = dat$ncdstype
    ),
    A = rep(list(1), length(effects_ncd)),
    effects = effects_ncd,
    tag = "ncd"
  )

  hiv_design <- as.data.frame(dat[, hiv_covariates, drop = FALSE])
  names(hiv_design) <- paste0(names(hiv_design), "_hiv")

  effects_hiv <- c(
    list(
      intercept_hiv = rep(1, n),
      region_id_hiv = dat$region_id,
      year_id_hiv = dat$year_id
    ),
    as.list(hiv_design)
  )

  stack_hiv <- INLA::inla.stack(
    data = list(
      response = cbind(NA_real_, dat$hiv_prev_p),
      log_expected = rep(0, n),
      country = dat$country,
      year = dat$year,
      ncdstype_diag = factor(NA, levels = levels(dat$ncdstype))
    ),
    A = rep(list(1), length(effects_hiv)),
    effects = effects_hiv,
    tag = "hiv"
  )

  stk <- INLA::inla.stack(stack_ncd, stack_hiv)

  ncd_fixed <- if (ncol(ncd_interactions) > 0) paste(colnames(ncd_interactions), collapse = " + ") else NULL
  hiv_fixed <- if (ncol(hiv_design) > 0) paste(colnames(hiv_design), collapse = " + ") else NULL

  fixed_terms <- c(
    "0 + intercept_ncd + ncdstype",
    ncd_fixed,
    "0 + intercept_hiv",
    hiv_fixed
  )

  formula_primary <- as.formula(paste0(
    "response ~ ", paste(fixed_terms[!vapply(fixed_terms, is.null, logical(1))], collapse = " + "), " + ",
    "f(region_id_ncd, model = 'bym2', graph = graph, hyper = bym2_hyper) + ",
    "f(year_id_ncd, model = 'rw1', hyper = rw1_prec) + ",
    "f(space_time_ncd, model = 'iid', hyper = iid_prec) + ",
    "f(region_id_hiv, model = 'bym2', graph = graph, hyper = bym2_hyper) + ",
    "f(year_id_hiv, model = 'rw1', hyper = rw1_prec)"
  ))
  environment(formula_primary) <- environment()

  stack_data <- INLA::inla.stack.data(stk)
  idx_ncd <- INLA::inla.stack.index(stk, tag = "ncd")$data

  exposure <- rep(NA_real_, nrow(stack_data$response))
  exposure[idx_ncd] <- exp(stack_data$log_expected[idx_ncd])

  result <- INLA::inla(
    formula_primary,
    family = c("poisson", "beta"),
    data = stack_data,
    E = exposure,
    control.family = list(
      list(link = "log"),
      list(link = "logit")
    ),
    control.predictor = list(A = INLA::inla.stack.A(stk), compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
    control.inla = list(strategy = "simplified.laplace")
  )

  list(res = result, stk = stk)
}

cat("\nFitting primary stacked-likelihood models...\n")
primary_fits <- lapply(names(model_tiers), function(model_name) {
  cat("  ", model_name, "\n", sep = "")
  build_primary_stacked_model(
    dat = dat,
    graph = graph,
    ncd_covariates = model_tiers[[model_name]],
    hiv_covariates = model_tiers[[model_name]]
  )
})
names(primary_fits) <- names(model_tiers)

saveRDS(primary_fits, file.path(model_dir, "primary_stacked_models.rds"))

primary_fit_table <- tibble::tibble(
  Model = names(primary_fits),
  DIC = vapply(primary_fits, function(x) x$res$dic$dic, numeric(1)),
  WAIC = vapply(primary_fits, function(x) x$res$waic$waic, numeric(1)),
  mean_logCPO = vapply(
    primary_fits,
    function(x) mean(log(pmax(x$res$cpo$cpo, .Machine$double.xmin)), na.rm = TRUE),
    numeric(1)
  )
) |>
  dplyr::mutate(
    delta_DIC = DIC - min(DIC, na.rm = TRUE),
    delta_WAIC = WAIC - min(WAIC, na.rm = TRUE)
  )

readr::write_csv(primary_fit_table, file.path(table_dir, "table_1_primary_model_fit.csv"))
print(primary_fit_table)

primary_M1 <- primary_fits$M1
primary_M4 <- primary_fits$M4
saveRDS(primary_M1, file.path(model_dir, "primary_M1.rds"))
saveRDS(primary_M4, file.path(model_dir, "primary_M4.rds"))

# Raw hyperparameter summaries are saved for traceability of Table 2.
primary_hyperpar <- dplyr::bind_rows(
  as.data.frame(primary_M1$res$summary.hyperpar) |>
    tibble::rownames_to_column("parameter") |>
    dplyr::mutate(Model = "M1"),
  as.data.frame(primary_M4$res$summary.hyperpar) |>
    tibble::rownames_to_column("parameter") |>
    dplyr::mutate(Model = "M4")
)
readr::write_csv(primary_hyperpar, file.path(table_dir, "table_2_primary_hyperparameters_raw.csv"))

# HIV fixed effects table for the primary model.
extract_hiv_fixed_effects <- function(model_object, model_label) {
  fixed <- as.data.frame(model_object$res$summary.fixed) |>
    tibble::rownames_to_column("term") |>
    dplyr::rename(
      mean = mean,
      lower = `0.025quant`,
      upper = `0.975quant`
    ) |>
    dplyr::filter(term == "intercept_hiv" | stringr::str_detect(term, "_hiv$")) |>
    dplyr::mutate(
      Model = model_label,
      Covariate = dplyr::case_when(
        term == "intercept_hiv" ~ "Intercept",
        term == "gdpercpta_z_hiv" ~ "GDP per capita, z",
        term == "avrgtemp_z_hiv" ~ "Mean annual temperature, z",
        term == "avrgprecp_z_hiv" ~ "Mean annual precipitation, z",
        term == "urbrte_z_hiv" ~ "Urbanisation rate, z",
        term == "hcexp_z_hiv" ~ "Health expenditure, z",
        TRUE ~ term
      ),
      OR = exp(mean),
      OR_L = exp(lower),
      OR_U = exp(upper)
    ) |>
    dplyr::select(Model, Covariate, mean, lower, upper, OR, OR_L, OR_U)

  fixed
}

hiv_fixed_effects <- dplyr::bind_rows(
  extract_hiv_fixed_effects(primary_M1, "M1"),
  extract_hiv_fixed_effects(primary_M4, "M4")
)
readr::write_csv(hiv_fixed_effects, file.path(table_dir, "table_4_hiv_fixed_effects.csv"))

# =============================================================================
# 5. Shared-component sensitivity model
# =============================================================================

safe_name <- function(x) {
  gsub("\\.", "_", make.names(as.character(x)))
}

build_shared_component_model <- function(dat, graph, ncd_covariates, hiv_covariates, model_label) {
  n <- nrow(dat)

  disease_levels <- levels(dat$ncdstype)
  disease_safe <- safe_name(disease_levels)

  disease_intercepts <- as.data.frame(
    sapply(seq_along(disease_levels), function(j) as.numeric(dat$ncdstype == disease_levels[j]))
  )
  names(disease_intercepts) <- paste0("ncd_int_", disease_safe)

  make_ncd_interactions <- function(df, covariates) {
    if (length(covariates) == 0) return(data.frame())

    output <- list()
    for (covariate in covariates) {
      for (j in seq_along(disease_levels)) {
        nm <- paste0("ncd_", safe_name(covariate), "__", disease_safe[j])
        output[[nm]] <- ifelse(df$ncdstype == disease_levels[j], df[[covariate]], 0)
      }
    }
    as.data.frame(output)
  }

  ncd_interactions <- make_ncd_interactions(dat, ncd_covariates)

  hiv_covariate_data <- as.data.frame(dat[, hiv_covariates, drop = FALSE])
  names(hiv_covariate_data) <- paste0("hiv_", safe_name(names(hiv_covariate_data)))

  fixed_cols <- c(
    "intercept_ncd", names(disease_intercepts), names(ncd_interactions),
    "intercept_hiv", names(hiv_covariate_data)
  )

  random_cols <- c(
    "region_shared", "region_shared_hiv",
    "year_shared", "year_shared_hiv",
    "space_time_ncd"
  )

  effects_ncd <- data.frame(
    intercept_ncd = rep(1, n),
    disease_intercepts,
    ncd_interactions,
    intercept_hiv = rep(0, n),
    hiv_covariate_data * 0,
    region_shared = dat$region_id,
    region_shared_hiv = NA_integer_,
    year_shared = dat$year_id,
    year_shared_hiv = NA_integer_,
    space_time_ncd = dat$space_time_id
  )

  effects_hiv <- data.frame(
    intercept_ncd = rep(0, n),
    disease_intercepts * 0,
    ncd_interactions * 0,
    intercept_hiv = rep(1, n),
    hiv_covariate_data,
    region_shared = NA_integer_,
    region_shared_hiv = dat$region_id,
    year_shared = NA_integer_,
    year_shared_hiv = dat$year_id,
    space_time_ncd = NA_integer_
  )

  for (nm in c(fixed_cols, random_cols)) {
    if (!nm %in% names(effects_ncd)) effects_ncd[[nm]] <- 0
    if (!nm %in% names(effects_hiv)) effects_hiv[[nm]] <- 0
  }

  effects_ncd <- effects_ncd[, c(fixed_cols, random_cols)]
  effects_hiv <- effects_hiv[, c(fixed_cols, random_cols)]

  stopifnot(nrow(effects_ncd) == n, nrow(effects_hiv) == n)
  stopifnot(identical(names(effects_ncd), names(effects_hiv)))

  stack_ncd <- INLA::inla.stack(
    data = list(
      response = cbind(dat$observed_deaths, NA_real_),
      log_expected = dat$log_expected,
      country = dat$country,
      year = dat$year,
      ncdstype_diag = dat$ncdstype
    ),
    A = list(1),
    effects = list(effects_ncd),
    tag = "ncd"
  )

  stack_hiv <- INLA::inla.stack(
    data = list(
      response = cbind(NA_real_, dat$hiv_prev_p),
      log_expected = rep(0, n),
      country = dat$country,
      year = dat$year,
      ncdstype_diag = factor(NA, levels = levels(dat$ncdstype))
    ),
    A = list(1),
    effects = list(effects_hiv),
    tag = "hiv"
  )

  stk <- INLA::inla.stack(stack_ncd, stack_hiv)

  copy_beta_prior <- list(
    beta = list(prior = "normal", param = c(0, 10), fixed = FALSE)
  )

  formula_shared <- as.formula(paste0(
    "response ~ 0 + ", paste(fixed_cols, collapse = " + "), " + ",
    "f(region_shared, model = 'bym2', graph = graph, hyper = bym2_hyper) + ",
    "f(region_shared_hiv, copy = 'region_shared', fixed = FALSE, hyper = copy_beta_prior) + ",
    "f(year_shared, model = 'rw1', hyper = rw1_prec) + ",
    "f(year_shared_hiv, copy = 'year_shared', fixed = FALSE, hyper = copy_beta_prior) + ",
    "f(space_time_ncd, model = 'iid', hyper = iid_prec)"
  ))
  environment(formula_shared) <- environment()

  stack_data <- INLA::inla.stack.data(stk)
  idx_ncd <- INLA::inla.stack.index(stk, tag = "ncd")$data
  idx_hiv <- INLA::inla.stack.index(stk, tag = "hiv")$data

  exposure <- rep(NA_real_, nrow(stack_data$response))
  exposure[idx_ncd] <- exp(stack_data$log_expected[idx_ncd])

  result <- INLA::inla(
    formula_shared,
    family = c("poisson", "beta"),
    data = stack_data,
    E = exposure,
    control.family = list(
      list(link = "log"),
      list(link = "logit")
    ),
    control.predictor = list(A = INLA::inla.stack.A(stk), compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
    control.inla = list(strategy = "simplified.laplace")
  )

  list(
    res = result,
    stk = stk,
    idx_ncd = idx_ncd,
    idx_hiv = idx_hiv,
    model_label = model_label,
    fixed_cols = fixed_cols
  )
}

cat("\nFitting shared-component sensitivity models...\n")
SC_M1 <- build_shared_component_model(dat, graph, M1_vars, M1_vars, "Shared-component M1")
SC_M4 <- build_shared_component_model(dat, graph, M4_vars, M4_vars, "Shared-component M4")

saveRDS(SC_M1, file.path(model_dir, "shared_component_M1.rds"))
saveRDS(SC_M4, file.path(model_dir, "shared_component_M4.rds"))

# =============================================================================
# 6. Shared-component summary tables
# =============================================================================

get_model_fit <- function(model_object, label) {
  tibble::tibble(
    Model = label,
    DIC = model_object$res$dic$dic,
    WAIC = model_object$res$waic$waic,
    mean_logCPO = mean(log(pmax(model_object$res$cpo$cpo, .Machine$double.xmin)), na.rm = TRUE)
  )
}

shared_model_fit <- dplyr::bind_rows(
  get_model_fit(SC_M1, "Shared-component M1"),
  get_model_fit(SC_M4, "Shared-component M4")
) |>
  dplyr::mutate(
    delta_DIC = DIC - min(DIC, na.rm = TRUE),
    delta_WAIC = WAIC - min(WAIC, na.rm = TRUE)
  )

readr::write_csv(shared_model_fit, file.path(table_dir, "shared_component_model_fit.csv"))

combined_model_fit <- dplyr::bind_rows(
  primary_fit_table |>
    dplyr::select(Model, DIC, WAIC, mean_logCPO) |>
    dplyr::mutate(Model = paste0("Primary stacked ", Model)),
  shared_model_fit |>
    dplyr::select(Model, DIC, WAIC, mean_logCPO)
) |>
  dplyr::mutate(
    delta_DIC = DIC - min(DIC, na.rm = TRUE),
    delta_WAIC = WAIC - min(WAIC, na.rm = TRUE)
  )

readr::write_csv(combined_model_fit, file.path(table_dir, "combined_model_fit.csv"))

extract_shared_loadings <- function(model_object, label) {
  hp <- as.data.frame(model_object$res$summary.hyperpar) |>
    tibble::rownames_to_column("parameter") |>
    dplyr::rename(
      mean = mean,
      lower = `0.025quant`,
      upper = `0.975quant`
    )

  hp |>
    dplyr::filter(stringr::str_detect(parameter, "Beta|beta|Copy|copy")) |>
    dplyr::mutate(
      component = dplyr::case_when(
        stringr::str_detect(parameter, "region_shared_hiv|region") ~ "HIV loading on shared BYM2 spatial field",
        stringr::str_detect(parameter, "year_shared_hiv|year") ~ "HIV loading on shared RW1 temporal field",
        TRUE ~ parameter
      ),
      Model = label
    ) |>
    dplyr::select(Model, component, parameter, mean, lower, upper)
}

shared_loadings <- dplyr::bind_rows(
  extract_shared_loadings(SC_M1, "Shared-component M1"),
  extract_shared_loadings(SC_M4, "Shared-component M4")
)

readr::write_csv(shared_loadings, file.path(table_dir, "shared_component_loadings.csv"))

extract_positive_probability <- function(model_object, label) {
  marginal_names <- names(model_object$res$marginals.hyperpar)
  beta_names <- marginal_names[stringr::str_detect(marginal_names, "Beta|beta|Copy|copy")]

  purrr::map_dfr(beta_names, function(parameter_name) {
    marginal <- model_object$res$marginals.hyperpar[[parameter_name]]
    tibble::tibble(
      Model = label,
      parameter = parameter_name,
      posterior_prob_positive = 1 - INLA::inla.pmarginal(0, marginal)
    )
  })
}

loading_probabilities <- dplyr::bind_rows(
  extract_positive_probability(SC_M1, "Shared-component M1"),
  extract_positive_probability(SC_M4, "Shared-component M4")
)

readr::write_csv(loading_probabilities, file.path(table_dir, "shared_loading_probability_positive.csv"))

make_table_3 <- function(shared_model_fit, shared_loadings, loading_probabilities) {
  loading_summary <- shared_loadings |>
    dplyr::filter(component %in% c(
      "HIV loading on shared BYM2 spatial field",
      "HIV loading on shared RW1 temporal field"
    )) |>
    dplyr::mutate(
      estimate = sprintf("%.3f [%.3f, %.3f]", mean, lower, upper),
      loading_type = dplyr::case_when(
        component == "HIV loading on shared BYM2 spatial field" ~ "spatial",
        component == "HIV loading on shared RW1 temporal field" ~ "temporal"
      )
    ) |>
    dplyr::select(Model, loading_type, estimate)

  probability_summary <- loading_probabilities |>
    dplyr::filter(stringr::str_detect(parameter, "region_shared_hiv|year_shared_hiv")) |>
    dplyr::mutate(
      loading_type = dplyr::case_when(
        stringr::str_detect(parameter, "region_shared_hiv") ~ "spatial",
        stringr::str_detect(parameter, "year_shared_hiv") ~ "temporal"
      )
    ) |>
    dplyr::select(Model, loading_type, posterior_prob_positive)

  loading_wide <- loading_summary |>
    tidyr::pivot_wider(names_from = loading_type, values_from = estimate, names_prefix = "loading_")

  probability_wide <- probability_summary |>
    tidyr::pivot_wider(
      names_from = loading_type,
      values_from = posterior_prob_positive,
      names_prefix = "Pr_"
    )

  shared_model_fit |>
    dplyr::select(Model, DIC, WAIC) |>
    dplyr::left_join(loading_wide, by = "Model") |>
    dplyr::left_join(probability_wide, by = "Model") |>
    dplyr::mutate(
      DIC = round(DIC, 0),
      WAIC = round(WAIC, 0),
      Pr_spatial = round(Pr_spatial, 3),
      Pr_temporal = round(Pr_temporal, 3)
    ) |>
    dplyr::rename(
      `HIV loading on shared BYM2 spatial field, mean [95% CrI]` = loading_spatial,
      `Pr(spatial loading > 0)` = Pr_spatial,
      `HIV loading on shared RW1 temporal field, mean [95% CrI]` = loading_temporal,
      `Pr(temporal loading > 0)` = Pr_temporal
    )
}

table_3_shared_component <- make_table_3(shared_model_fit, shared_loadings, loading_probabilities)
readr::write_csv(table_3_shared_component, file.path(table_dir, "table_3_shared_component_sensitivity.csv"))
print(table_3_shared_component)

# =============================================================================
# 7. Shared spatial component map: manuscript Figure 2
# =============================================================================

extract_shared_spatial_sf <- function(model_object, shp) {
  if (!"region_shared" %in% names(model_object$res$summary.random)) {
    stop("The model does not contain a 'region_shared' random effect.", call. = FALSE)
  }

  spatial_df <- as.data.frame(model_object$res$summary.random$region_shared) |>
    dplyr::rename(
      region_id = ID,
      mean = mean,
      lower = `0.025quant`,
      upper = `0.975quant`
    ) |>
    dplyr::mutate(region_id = as.integer(region_id))

  shp |>
    dplyr::mutate(region_id = as.integer(region_id)) |>
    dplyr::left_join(spatial_df, by = "region_id")
}

sf_SC_M1 <- extract_shared_spatial_sf(SC_M1, shp) |>
  dplyr::mutate(Model = "(A) Shared-component M1")

sf_SC_M4 <- extract_shared_spatial_sf(SC_M4, shp) |>
  dplyr::mutate(Model = "(B) Shared-component M4")

sf_shared_spatial <- rbind(
  sf_SC_M1 |> dplyr::select(country, region_id, Model, mean, lower, upper, geometry),
  sf_SC_M4 |> dplyr::select(country, region_id, Model, mean, lower, upper, geometry)
)

sf_shared_spatial$Model <- factor(
  sf_shared_spatial$Model,
  levels = c("(A) Shared-component M1", "(B) Shared-component M4")
)

spatial_breaks <- c(-0.8, -0.4, -0.2, -0.1, 0, 0.1, 0.2, 0.4)

tmap::tmap_mode("plot")

figure_2_shared_spatial <- tmap::tm_shape(sf_shared_spatial) +
  tmap::tm_polygons(
    fill = "mean",
    fill.scale = tmap::tm_scale_intervals(
      values = "-brewer.rd_yl_bu",
      breaks = spatial_breaks,
      midpoint = 0
    ),
    fill.legend = tmap::tm_legend(
      title = "Shared spatial effect",
      title.fontface = "bold",
      orientation = "landscape",
      title.position = "top",
      title.align = "center",
      width = 25.6,
      height = 2.8,
      text.size = 1.1,
      frame = FALSE
    )
  ) +
  tmap::tm_borders(col = "black", lwd = 0.7) +
  tmap::tm_facets(by = "Model", ncol = 2) +
  tmap::tm_compass(type = "rose", position = c("right", "top"), size = 2.5) +
  tmap::tm_title(
    "Shared BYM2 spatial component linking NCD mortality and HIV prevalence",
    size = 1.1,
    fontface = "bold"
  ) +
  tmap::tm_layout(
    legend.outside = TRUE,
    legend.outside.position = "bottom",
    legend.title.size = 0.9,
    legend.text.size = 0.7,
    panel.label.size = 1.1,
    outer.margins = c(0.02, 0.02, 0.02, 0.02),
    inner.margins = c(0.12, 0.06, 0.06, 0.06)
  )

figure_2_shared_spatial

tmap::tmap_save(
  figure_2_shared_spatial,
  file.path(figure_dir, "figure_2_shared_spatial_component.png"),
  dpi = 1000,
  width = 7.5,
  height = 5.5,
  units = "in"
)

# =============================================================================
# 8. Shared temporal component: manuscript Figure 3
# =============================================================================

extract_shared_temporal <- function(model_object, label, years = 2000:2019) {
  if (!"year_shared" %in% names(model_object$res$summary.random)) {
    stop("The model does not contain a 'year_shared' random effect.", call. = FALSE)
  }

  temporal_df <- as.data.frame(model_object$res$summary.random$year_shared)

  tibble::tibble(
    year = years[seq_len(nrow(temporal_df))],
    mean = temporal_df$mean,
    lower = temporal_df$`0.025quant`,
    upper = temporal_df$`0.975quant`,
    Model = label
  )
}

shared_temporal <- dplyr::bind_rows(
  extract_shared_temporal(SC_M1, "Shared-component M1"),
  extract_shared_temporal(SC_M4, "Shared-component M4")
)

readr::write_csv(shared_temporal, file.path(table_dir, "shared_temporal_effects.csv"))

figure_3_shared_temporal <- ggplot2::ggplot(
  shared_temporal,
  ggplot2::aes(x = year, y = mean, colour = Model, fill = Model)
) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA) +
  ggplot2::geom_line(linewidth = 1.05) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, colour = "grey40") +
  ggplot2::labs(
    title = "Shared RW1 temporal component linking NCD mortality and HIV prevalence",
    x = "Year",
    y = "Posterior mean shared temporal effect"
  ) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.08))) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5, margin = ggplot2::margin(b = 12)),
    axis.title.x = ggplot2::element_text(size = 12.5, margin = ggplot2::margin(t = 12)),
    axis.title.y = ggplot2::element_text(size = 12.5, margin = ggplot2::margin(r = 12)),
    axis.text = ggplot2::element_text(size = 11, colour = "black"),
    panel.grid.major = ggplot2::element_line(colour = "grey85", linewidth = 0.4),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 11),
    plot.margin = ggplot2::margin(t = 10, r = 15, b = 10, l = 15)
  )

figure_3_shared_temporal

ggplot2::ggsave(
  file.path(figure_dir, "figure_3_shared_temporal_component.png"),
  figure_3_shared_temporal,
  width = 7.6,
  height = 6,
  dpi = 800
)

# =============================================================================
# 9. Residual co-clustering after shared-component adjustment: Figure 4
# =============================================================================

compute_residual_correlation_shared <- function(model_object, dat, label) {
  idx_ncd <- INLA::inla.stack.index(model_object$stk, tag = "ncd")$data
  idx_hiv <- INLA::inla.stack.index(model_object$stk, tag = "hiv")$data
  fitted_values <- model_object$res$summary.fitted.values$mean

  dat |>
    dplyr::mutate(
      resid_ncd = observed_deaths - fitted_values[idx_ncd],
      resid_hiv = hiv_prev_p - fitted_values[idx_hiv]
    ) |>
    dplyr::group_by(region_id, country) |>
    dplyr::summarise(
      residual_correlation = stats::cor(resid_ncd, resid_hiv, use = "complete.obs"),
      .groups = "drop"
    ) |>
    dplyr::mutate(Model = label)
}

residual_correlation_shared <- dplyr::bind_rows(
  compute_residual_correlation_shared(SC_M1, dat, "(A) Shared-component M1"),
  compute_residual_correlation_shared(SC_M4, dat, "(B) Shared-component M4")
)

readr::write_csv(residual_correlation_shared, file.path(table_dir, "residual_correlation_shared_component.csv"))

shp_residual <- shp |>
  dplyr::select(-dplyr::any_of("region_id")) |>
  dplyr::left_join(dat |> dplyr::distinct(country, region_id), by = "country") |>
  dplyr::left_join(residual_correlation_shared, by = c("country", "region_id"))

if (any(is.na(shp_residual$residual_correlation))) {
  stop("Residual correlations did not join completely to the spatial data.", call. = FALSE)
}

residual_limit <- max(abs(shp_residual$residual_correlation), na.rm = TRUE)
residual_breaks <- seq(
  -ceiling(residual_limit * 10) / 10,
  ceiling(residual_limit * 10) / 10,
  by = 0.2
)

figure_4_residual_map <- tmap::tm_shape(shp_residual) +
  tmap::tm_polygons(
    fill = "residual_correlation",
    fill.scale = tmap::tm_scale_intervals(
      values = "-brewer.rd_bu",
      breaks = residual_breaks,
      midpoint = 0
    ),
    fill.legend = tmap::tm_legend(
      title = "Residual correlation",
      title.fontface = "bold",
      orientation = "landscape",
      title.position = "top",
      title.align = "center",
      width = 10.6,
      height = 2.8,
      text.size = 1,
      title.padding = c(0.01, 0.01, 0.05, 0.01),
      frame = FALSE
    )
  ) +
  tmap::tm_borders(col = "black", lwd = 0.7) +
  tmap::tm_facets(by = "Model", ncol = 2) +
  tmap::tm_compass(type = "rose", position = c("right", "top"), size = 2.5) +
  tmap::tm_title(
    "Residual HIV-NCD co-clustering after shared-component adjustment",
    size = 1.1,
    fontface = "bold"
  ) +
  tmap::tm_layout(
    legend.outside = TRUE,
    legend.outside.position = "bottom",
    legend.title.size = 0.9,
    legend.text.size = 0.7,
    panel.label.size = 1.1,
    outer.margins = c(0.02, 0.02, 0.02, 0.02),
    inner.margins = c(0.12, 0.06, 0.06, 0.06)
  )

figure_4_residual_map

tmap::tmap_save(
  figure_4_residual_map,
  filename = file.path(figure_dir, "figure_4_residual_coclustering_shared_component.png"),
  dpi = 1000,
  width = 7.5,
  height = 5.5,
  units = "in"
)

# =============================================================================
# 10. HIV calibration under shared-component model: Figure 5
# =============================================================================

make_hiv_observed_fitted <- function(model_object, dat, label) {
  idx_hiv <- INLA::inla.stack.index(model_object$stk, tag = "hiv")$data
  fitted_hiv <- model_object$res$summary.fitted.values$mean[idx_hiv]

  tibble::tibble(
    country = dat$country,
    year = dat$year,
    observed = dat$hiv_prev,
    fitted = fitted_hiv * 100,
    Model = label
  )
}

hiv_observed_fitted <- dplyr::bind_rows(
  make_hiv_observed_fitted(SC_M1, dat, "Shared-component M1"),
  make_hiv_observed_fitted(SC_M4, dat, "Shared-component M4")
)

readr::write_csv(hiv_observed_fitted, file.path(table_dir, "hiv_observed_fitted_shared_component.csv"))

hiv_panel_stats <- hiv_observed_fitted |>
  dplyr::group_by(Model) |>
  dplyr::summarise(
    R2 = stats::cor(observed, fitted, use = "complete.obs")^2,
    RMSE = sqrt(mean((observed - fitted)^2, na.rm = TRUE)),
    x = min(fitted, na.rm = TRUE) + 0.05 * diff(range(fitted, na.rm = TRUE)),
    y = max(observed, na.rm = TRUE) - 0.08 * diff(range(observed, na.rm = TRUE)),
    label = sprintf("R² = %.2f\nRMSE = %.2f", R2, RMSE),
    .groups = "drop"
  )

figure_5_hiv_calibration <- ggplot2::ggplot(
  hiv_observed_fitted,
  ggplot2::aes(x = fitted, y = observed, colour = country)
) +
  ggplot2::geom_point(alpha = 0.75, size = 2.0) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.6, colour = "grey40") +
  ggplot2::geom_text(
    data = hiv_panel_stats,
    ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.5
  ) +
  ggplot2::facet_wrap(~ Model) +
  ggplot2::labs(
    title = "Observed versus fitted HIV prevalence under the shared-component model",
    x = "Fitted HIV prevalence (%)",
    y = "Observed HIV prevalence (%)",
    colour = "Country"
  ) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.03))) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.05))) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 13.6, hjust = 0.5, margin = ggplot2::margin(b = 12)),
    axis.title.x = ggplot2::element_text(size = 12.5, margin = ggplot2::margin(t = 12)),
    axis.title.y = ggplot2::element_text(size = 12.5, margin = ggplot2::margin(r = 12)),
    axis.text = ggplot2::element_text(size = 11, colour = "black"),
    strip.text = ggplot2::element_text(size = 12.5, face = "bold"),
    panel.grid.major = ggplot2::element_line(colour = "grey85", linewidth = 0.4),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 12),
    legend.text = ggplot2::element_text(size = 10),
    plot.margin = ggplot2::margin(t = 10, r = 15, b = 10, l = 15)
  )

figure_5_hiv_calibration

ggplot2::ggsave(
  file.path(figure_dir, "figure_5_hiv_observed_fitted_shared_component.png"),
  figure_5_hiv_calibration,
  width = 7.5,
  height = 6,
  dpi = 800
)

# =============================================================================
# 11. NCD calibration under shared-component model: Figure 6
# =============================================================================

make_ncd_observed_fitted <- function(model_object, dat, label) {
  idx_ncd <- INLA::inla.stack.index(model_object$stk, tag = "ncd")$data

  # For Poisson models fitted with the E argument, fitted values are relative
  # risks/rates. Fitted counts are obtained by multiplying by expected deaths.
  fitted_relative_risk <- model_object$res$summary.fitted.values$mean[idx_ncd]
  expected_deaths <- exp(dat$log_expected)
  fitted_counts <- fitted_relative_risk * expected_deaths

  tibble::tibble(
    country = dat$country,
    year = dat$year,
    disease = dat$ncdstype,
    observed = dat$observed_deaths,
    expected = expected_deaths,
    fitted_relative_risk = fitted_relative_risk,
    fitted = fitted_counts,
    Model = label
  )
}

ncd_observed_fitted <- dplyr::bind_rows(
  make_ncd_observed_fitted(SC_M1, dat, "Shared-component M1"),
  make_ncd_observed_fitted(SC_M4, dat, "Shared-component M4")
)

readr::write_csv(ncd_observed_fitted, file.path(table_dir, "ncd_observed_fitted_shared_component.csv"))

ncd_fit_stats <- ncd_observed_fitted |>
  dplyr::group_by(Model, disease) |>
  dplyr::summarise(
    n = dplyr::n(),
    correlation = stats::cor(observed, fitted, use = "complete.obs"),
    MAE = mean(abs(observed - fitted), na.rm = TRUE),
    RMSE = sqrt(mean((observed - fitted)^2, na.rm = TRUE)),
    mean_observed = mean(observed, na.rm = TRUE),
    mean_fitted = mean(fitted, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(ncd_fit_stats, file.path(table_dir, "ncd_fit_diagnostics_shared_component.csv"))

disease_labels <- c(
  cvds = "Cardiovascular\ndiseases",
  diabmtus = "Diabetes\nmellitus",
  malineopms = "Malignant\nneoplasms",
  rsptoryds = "Respiratory\ndiseases"
)

log_breaks <- c(1e3, 1e4, 1e5)

figure_6_ncd_calibration <- ggplot2::ggplot(
  ncd_observed_fitted,
  ggplot2::aes(x = fitted + 1, y = observed + 1, colour = country)
) +
  ggplot2::geom_point(alpha = 0.7, size = 1.8) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.5, colour = "grey40") +
  ggplot2::scale_x_log10(breaks = log_breaks, labels = scales::label_comma()) +
  ggplot2::scale_y_log10(breaks = log_breaks, labels = scales::label_comma()) +
  ggplot2::facet_grid(Model ~ disease, labeller = ggplot2::labeller(disease = disease_labels)) +
  ggplot2::labs(
    title = "Observed versus fitted NCD deaths under the shared-component model",
    x = "Fitted NCD deaths (log10 scale)",
    y = "Observed NCD deaths (log10 scale)",
    colour = "Country"
  ) +
  ggplot2::theme_minimal(base_size = 12.5) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 13.6, hjust = 0.5, margin = ggplot2::margin(b = 12)),
    axis.title.x = ggplot2::element_text(size = 12, margin = ggplot2::margin(t = 14)),
    axis.title.y = ggplot2::element_text(size = 12, margin = ggplot2::margin(r = 14)),
    axis.text = ggplot2::element_text(size = 10.5, colour = "black"),
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    strip.text.x = ggplot2::element_text(face = "bold", size = 9.5),
    strip.text.y = ggplot2::element_text(face = "bold", size = 11.5),
    panel.grid.major = ggplot2::element_line(colour = "grey85", linewidth = 0.35),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 11.5),
    legend.text = ggplot2::element_text(size = 9.5),
    panel.spacing = grid::unit(1.1, "lines"),
    plot.margin = ggplot2::margin(t = 10, r = 16, b = 10, l = 16)
  )

figure_6_ncd_calibration

ggplot2::ggsave(
  filename = file.path(figure_dir, "figure_6_ncd_observed_fitted_shared_component.png"),
  plot = figure_6_ncd_calibration,
  width = 7.5,
  height = 6.5,
  dpi = 1000
)

cat("\nAnalysis complete. Outputs saved in: ", out_dir, "\n", sep = "")
