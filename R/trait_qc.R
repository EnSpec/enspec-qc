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

# ---- 5. Assay replicate precision (relative error) ------------------------

#' Flag assay results whose own replicate scatter is too large relative to
#' the value -- i.e. the measurement disagreed with itself.
#'
#' Distinct from every other check here: it needs no comparison to other
#' samples, other groups or any literature range. An assay run in triplicate
#' that reports a mean and an SD carries its own precision estimate, and a
#' coefficient of variation above what the method can justify means that
#' result is unreliable regardless of how plausible the value looks.
#'
#' @param value_col the assay result (a mean over replicates)
#' @param error_col the matching dispersion (SD across those replicates)
#' @param max_cv flag when error_col / value_col exceeds this. Set it from
#'   what the assay can actually achieve, per project.
qc_relative_error <- function(df, value_col, error_col, id_cols, max_cv = 0.10) {
  df %>%
    filter(!is.na(.data[[value_col]]), !is.na(.data[[error_col]]),
           .data[[value_col]] > 0) %>%
    mutate(cv = .data[[error_col]] / .data[[value_col]]) %>%
    filter(cv > max_cv) %>%
    select(all_of(id_cols), all_of(value_col), all_of(error_col), cv) %>%
    mutate(trait_flagged = value_col,
           max_cv = max_cv,
           check = "relative_error") %>%
    arrange(desc(cv))
}

# ---- 6. QC blank drift across sequential processing steps -----------------

#' Flag a QC blank that GAINS weight across a sequence of processing steps.
#'
#' For a sequential gravimetric assay (ANKOM fiber being the motivating
#' case: NDF then ADF then ADL), an empty blank bag carried through the same
#' washes should only ever lose weight. If it gains, the bag leaked and took
#' on material -- which means it was shedding sample material too, and
#' everything sharing that batch is suspect. This is a BATCH-level signal,
#' not a per-sample one: the consequence of a leaking blank falls on every
#' sample processed alongside it, so the caller is responsible for
#' propagating the flag to the batch.
#'
#' Reports both comparisons, because they answer different questions:
#'   - change from tare: did the blank end up heavier than it started
#'   - change from the previous step: which wash was it that added weight
#'
#' @param tare_col the blank's initial empty weight
#' @param step_cols weights after each step, IN PROCESSING ORDER
#' @param increase_threshold fractional gain treated as bad data (0.01 = 1%,
#'   a threshold in common use). Any gain at all is reported regardless, since
#'   a blank should not gain weight even slightly.
qc_blank_drift <- function(df, id_cols, tare_col, step_cols,
                            increase_threshold = 0.01) {
  purrr::map_dfr(seq_len(nrow(df)), function(i) {
    tare <- df[[tare_col]][i]
    if (is.na(tare) || tare <= 0) return(NULL)

    weights <- vapply(step_cols, function(cl) as.numeric(df[[cl]][i]), numeric(1))
    prev <- c(tare, weights[-length(weights)])

    tibble::tibble(
      df[i, id_cols, drop = FALSE],
      step = step_cols,
      step_order = seq_along(step_cols),
      weight = weights,
      tare = tare,
      pct_change_from_tare = 100 * (weights / tare - 1),
      pct_change_from_previous = 100 * (weights / prev - 1)
    )
  }) %>%
    mutate(
      gained_weight = pct_change_from_tare > 0 | pct_change_from_previous > 0,
      exceeds_threshold = pmax(pct_change_from_tare, pct_change_from_previous, na.rm = TRUE) >
        100 * increase_threshold,
      increase_threshold_pct = 100 * increase_threshold,
      check = "blank_drift"
    ) %>%
    filter(gained_weight) %>%
    arrange(desc(pmax(pct_change_from_tare, pct_change_from_previous, na.rm = TRUE)))
}

# ---- 7. Sequential-assay cascade -------------------------------------------

#' Propagate a failure forward through a sequential assay.
#'
#' Some assays produce their traits in a fixed order, each stage operating on
#' the residue of the last -- ANKOM fiber (NDF, then ADF on the NDF residue,
#' then ADL on the ADF residue) being the motivating case. When an early
#' stage fails, every trait derived after it is suspect, whether or not its
#' own value looks reasonable. Judging each trait in isolation misses that
#' entirely (Henry, 2026-09-23).
#'
#' This is deliberately precautionary: it flags on position in the chain, not
#' on evidence that the downstream value is itself wrong. Depending on how a
#' given assay's arithmetic actually wires up, a downstream trait may be
#' untouched by the upstream failure -- so treat these as "look at this",
#' and check the assay's formulas before removing anything on this basis.
#'
#' @param failed_flags flag table for the traits already known to have failed
#' @param chain trait names IN PROCESSING ORDER
#' @param key_cols columns identifying the unit the chain runs on (e.g.
#'   sample and treatment) -- the cascade stays within one unit
#' @param trait_col column in failed_flags naming the failed trait
qc_sequential_cascade <- function(failed_flags, chain, key_cols,
                                   trait_col = "trait_flagged") {
  if (nrow(failed_flags) == 0) return(failed_flags[0, ] %>% mutate(check = character(0)))

  in_chain <- failed_flags[failed_flags[[trait_col]] %in% chain, , drop = FALSE]
  if (nrow(in_chain) == 0) return(in_chain %>% mutate(check = character(0)))

  purrr::pmap_dfr(in_chain[, c(key_cols, trait_col)], function(...) {
    row <- list(...)
    failed_trait <- row[[trait_col]]
    pos <- match(failed_trait, chain)
    downstream <- chain[seq_len(length(chain)) > pos]
    if (length(downstream) == 0) return(NULL)

    tibble::tibble(
      !!!row[key_cols],
      trait_flagged = downstream,
      upstream_failure = failed_trait,
      position_in_chain = match(downstream, chain),
      reason = sprintf("downstream of failed '%s' in a sequential assay -- suspect by position, verify against the assay's own arithmetic before acting",
                       failed_trait),
      check = "sequential_cascade"
    )
  }) %>%
    distinct()
}

# ---- 8. Trait-trait covariation flag ---------------------------------------

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
