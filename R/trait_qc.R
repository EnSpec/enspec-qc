##################################################
#
# trait_qc.R
#
# Shared, project-agnostic trait QC functions.
#
# Design (decided for the Mill Test 2024 EcoSIS deposit, 2026-09-20,
# generalized here as the first cross-project implementation):
#   1. Hard-bounds sweep for impossible/mis-entered values -> qc_hard_bounds()
#   2. Bad single-measure / ratio-component checks (leave-one-out median,
#      catches decimal shifts and instrument mis-reads within a sampling
#      unit) -> qc_bad_single_measure(), qc_bad_ratio_component()
#   3. Magnitude/decimal-shift check against a group median -> qc_magnitude_shift()
#   4. Distribution-based flag (robust z, within group) -> qc_distribution_outliers()
#   5. Trait-trait covariation flag (robust SMA residuals) -> qc_covariation_outliers()
#
# Every function takes the "sampling unit" as a caller-supplied grouping
# variable (or variables, via a single group_var column you build upstream,
# e.g. paste(Species_Code, DryTreat)) -- what counts as "variability" is
# project-specific (taxon, treatment, treatment x taxon, ...); the
# comparison logic here is not.
#
# All functions FLAG by default (return the rows that tripped the check,
# they do not modify or drop anything). Whether a flag becomes a removal is
# a per-project policy decision, applied downstream with qc_apply_removals()
# in qc_report.R -- keeps this file itself free of any removal side effects.
#
##################################################

library(dplyr)
library(rlang)

#' Leave-one-out group median, aligned to x's original order
.loo_median <- function(x) {
  vapply(seq_along(x), function(i) median(x[-i], na.rm = TRUE), numeric(1))
}

# ---- 1. Hard-bounds sweep --------------------------------------------------

#' Flag values outside a physically-possible or literature-supported range.
#'
#' @param df data frame
#' @param bounds a data frame/tibble with columns: trait, lower, upper,
#'   source (citation or "physical" for e.g. percent 0-100), reason (free text)
#' @param id_cols character vector of columns to carry through for provenance
#' @return long-format tibble, one row per (id, trait) flag, with the
#'   offending value and the bound(s) it violated
qc_hard_bounds <- function(df, bounds, id_cols) {
  stopifnot(all(c("trait", "lower", "upper") %in% names(bounds)))

  purrr::pmap_dfr(bounds, function(trait, lower, upper, source = NA, reason = NA, ...) {
    if (!trait %in% names(df)) {
      warning(sprintf("qc_hard_bounds: trait '%s' not found in df, skipping", trait))
      return(NULL)
    }
    if (is.na(lower) && is.na(upper)) return(NULL)  # bounds not yet filled in

    vals <- df[[trait]]
    bad <- !is.na(vals) & ((!is.na(lower) & vals < lower) | (!is.na(upper) & vals > upper))
    if (!any(bad)) return(NULL)

    df[bad, id_cols, drop = FALSE] %>%
      mutate(
        trait_flagged = trait,
        value = vals[bad],
        lower_bound = lower,
        upper_bound = upper,
        bound_source = source,
        reason = reason,
        check = "hard_bounds"
      )
  })
}

# ---- 2. Bad single-measure / ratio-component checks ------------------------

#' Flag a single measured value that's off by a large factor from the
#' leave-one-out median of its own sampling unit (e.g. a replicate weight
#' that's 10x the others in its group -- a likely decimal-shift or
#' balance-entry error, not necessarily a "real" outlier).
#'
#' @param group_var column name (string) defining the sampling unit
#'   (e.g. "sample_id", or a combined "species_x_treatment" column you built)
#' @param ratio_threshold flag if value / loo_median > threshold or < 1/threshold
qc_bad_single_measure <- function(df, group_var, value_col, id_cols,
                                   ratio_threshold = 2, min_n = 2) {
  df %>%
    filter(!is.na(.data[[value_col]])) %>%
    group_by(.data[[group_var]]) %>%
    filter(n() >= min_n) %>%
    mutate(
      group_loo_median = .loo_median(.data[[value_col]]),
      value_ratio = .data[[value_col]] / group_loo_median,
      flag_bad_measure = value_ratio > ratio_threshold | value_ratio < (1 / ratio_threshold)
    ) %>%
    ungroup() %>%
    filter(flag_bad_measure) %>%
    select(all_of(id_cols), all_of(group_var), all_of(value_col),
           group_loo_median, value_ratio) %>%
    mutate(trait_flagged = value_col, check = "bad_single_measure") %>%
    arrange(desc(value_ratio))
}

#' Same idea, but for a trait computed as a ratio/difference of two raw
#' components (e.g. LMA = dry_wgt / area, %oil = oil_mass / sample_wgt).
#' Flags on the ratio trait itself, then attributes the likely bad
#' component by comparing each component's own leave-one-out deviation --
#' whichever component moved more is "suspected_bad".
qc_bad_ratio_component <- function(df, group_var, ratio_col, comp1_col, comp2_col,
                                    id_cols, ratio_threshold = 2, min_n = 2) {
  df %>%
    filter(!is.na(.data[[ratio_col]])) %>%
    group_by(.data[[group_var]]) %>%
    filter(n() >= min_n) %>%
    mutate(
      ratio_loo_med = .loo_median(.data[[ratio_col]]),
      comp1_loo_med = .loo_median(.data[[comp1_col]]),
      comp2_loo_med = .loo_median(.data[[comp2_col]]),

      trait_ratio = .data[[ratio_col]] / ratio_loo_med,
      comp1_dev = abs(log2(.data[[comp1_col]] / comp1_loo_med)),
      comp2_dev = abs(log2(.data[[comp2_col]] / comp2_loo_med)),

      flag_bad = trait_ratio > ratio_threshold | trait_ratio < (1 / ratio_threshold),
      suspected_bad = case_when(
        !flag_bad ~ NA_character_,
        comp1_dev > comp2_dev ~ comp1_col,
        comp2_dev > comp1_dev ~ comp2_col,
        TRUE ~ "ambiguous"
      )
    ) %>%
    ungroup() %>%
    filter(flag_bad) %>%
    select(all_of(id_cols), all_of(group_var),
           all_of(ratio_col), trait_ratio,
           all_of(comp1_col), comp1_dev,
           all_of(comp2_col), comp2_dev,
           suspected_bad) %>%
    mutate(trait_flagged = ratio_col, check = "bad_ratio_component") %>%
    arrange(desc(trait_ratio))
}

# ---- 3. Magnitude / decimal-shift check ------------------------------------

#' Flag values far (in log10 space) from their group's median -- catches
#' decimal-shift-style entry errors at a coarser grain than
#' qc_bad_single_measure (bigger groups, e.g. species or species x treatment,
#' not just replicate sets within one sample).
qc_magnitude_shift <- function(df, group_var, value_col, id_cols,
                                min_n = 5, log10_threshold = 1) {
  df %>%
    filter(!is.na(.data[[value_col]]), .data[[value_col]] > 0) %>%
    group_by(.data[[group_var]]) %>%
    filter(n() >= min_n) %>%
    mutate(
      group_median = median(.data[[value_col]], na.rm = TRUE),
      log10_dev = abs(log10(.data[[value_col]] / group_median)),
      flag_magnitude = log10_dev > log10_threshold
    ) %>%
    ungroup() %>%
    filter(flag_magnitude) %>%
    select(all_of(id_cols), all_of(group_var), all_of(value_col),
           group_median, log10_dev) %>%
    mutate(trait_flagged = value_col, check = "magnitude_shift") %>%
    arrange(desc(log10_dev))
}

# ---- 4. Distribution-based flag (robust, within group) ---------------------

#' Robust (MAD-based) z-score flag within a group. Flags, never removes --
#' distribution outliers can be real biological variation, not entry errors.
qc_distribution_outliers <- function(df, group_var, value_col, id_cols,
                                      min_n = 5, z_threshold = 3.5) {
  df %>%
    filter(!is.na(.data[[value_col]])) %>%
    group_by(.data[[group_var]]) %>%
    filter(n() >= min_n) %>%
    mutate(
      group_median = median(.data[[value_col]], na.rm = TRUE),
      group_mad = mad(.data[[value_col]], na.rm = TRUE),
      # 1.4826 * MAD ~= SD under normality; guard the degenerate all-equal case
      robust_z = if_else(group_mad > 0,
                         (.data[[value_col]] - group_median) / (1.4826 * group_mad),
                         0),
      flag_distribution = abs(robust_z) > z_threshold
    ) %>%
    ungroup() %>%
    filter(flag_distribution) %>%
    select(all_of(id_cols), all_of(group_var), all_of(value_col),
           group_median, robust_z) %>%
    mutate(trait_flagged = value_col, check = "distribution_outlier") %>%
    arrange(desc(abs(robust_z)))
}

# ---- 5. Trait-trait covariation flag ---------------------------------------

#' Flag points far from a robust (Huber M-estimation) SMA fit between two
#' traits expected to strongly covary (e.g. N vs LMA in the leaf economics
#' spectrum). Fit on log10-log10 scale by default -- set log_scale = FALSE
#' for traits where a linear relationship is more appropriate.
#' If group_var is supplied, the fit + residual z-scoring is done separately
#' within each group (use this when the covariation relationship itself is
#' expected to differ by group, e.g. by growth form); otherwise one global
#' fit is used (matches the GCFR precedent for N-LMA).
qc_covariation_outliers <- function(df, trait1, trait2, id_cols,
                                    group_var = NULL, log_scale = TRUE,
                                    z_threshold = 3, min_n = 8) {
  if (!requireNamespace("smatr", quietly = TRUE)) {
    stop("qc_covariation_outliers() requires the 'smatr' package")
  }

  fit_one <- function(sub) {
    sub <- sub %>% filter(!is.na(.data[[trait1]]), !is.na(.data[[trait2]]))
    if (log_scale) sub <- sub %>% filter(.data[[trait1]] > 0, .data[[trait2]] > 0)
    if (nrow(sub) < min_n) return(NULL)

    x <- if (log_scale) log10(sub[[trait2]]) else sub[[trait2]]
    y <- if (log_scale) log10(sub[[trait1]]) else sub[[trait1]]

    fit <- smatr::sma(y ~ x, robust = TRUE)
    resid_z <- as.numeric(scale(residuals(fit)))

    sub %>%
      mutate(sma_resid_z = resid_z,
             flag_covariation = abs(sma_resid_z) > z_threshold) %>%
      filter(flag_covariation)
  }

  result <- if (is.null(group_var)) {
    fit_one(df)
  } else {
    df %>% group_by(.data[[group_var]]) %>% group_modify(~ fit_one(.x) %||% .x[0, ]) %>% ungroup()
  }

  if (is.null(result) || nrow(result) == 0) {
    return(df[0, c(id_cols, trait1, trait2)] %>%
             mutate(sma_resid_z = numeric(0), trait_flagged = character(0), check = character(0)))
  }

  result %>%
    select(all_of(id_cols), any_of(group_var), all_of(trait1), all_of(trait2), sma_resid_z) %>%
    mutate(trait_flagged = paste(trait1, "vs", trait2), check = "covariation_outlier") %>%
    arrange(desc(abs(sma_resid_z)))
}
