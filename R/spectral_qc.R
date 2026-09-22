##################################################
#
# spectral_qc.R
#
# Shared, project-agnostic spectral QC functions for reflectance spectra
# measured with a contact probe / press apparatus on prepared (e.g. dried,
# milled) plant material -- not a field/canopy spectroscopy QC pipeline.
#
# Unlike trait QC, there's no existing internal precedent for this (first
# built for the Mill Test 2024 deposit, 2026-09). The measurement is fairly
# constrained by the method, so outliers should be rare in a controlled
# system -- keep thresholds conservative, and treat these as flags, not
# removals, until there's enough cross-project experience to say otherwise.
#
# Known limitation (Henry, 2026-09-20): the most useful check would be a
# moisture-index-style flag for incomplete drying, but that needs a
# validation dataset this project doesn't have -- not attempted here.
# What IS attempted: catching bad reference-puck reads / the probe not
# sitting fully on the sample, via two conservative proxies:
#   1. Flat-and-extreme spectra (a measurement stuck near 0 or near full
#      reflectance, with no real spectral shape) -> qc_spectral_flat()
#   2. Physically-impossible reflectance values -> qc_spectral_out_of_range()
# Plus one within-group check:
#   3. Spectra that are unusually dissimilar in shape from their own
#      group's typical spectrum -> qc_spectral_group_outliers()
#
##################################################

library(dplyr)

# ---- 1. Flat / extreme (bad reference or probe-off-sample) proxy ----------

#' Flag spectra that are both (a) nearly flat across wavelength (low
#' variance = "looks like a straight line") and (b) sitting at an extreme
#' reflectance level (near 0% or near 100%). A real leaf/tissue spectrum has
#' shape (absorption features); a flat line pinned near 0 or 1 across the
#' full range is much more consistent with the probe not making contact or
#' a bad white-reference read than with real material.
#'
#' @param wave_cols character vector of wavelength column names
#' @param flat_range_threshold flag as "flat" if (max-min) across wave_cols
#'   is below this (on the same 0-1 reflectance scale)
#' @param low_thresh / high_thresh median reflectance bounds defining "extreme"
qc_spectral_flat <- function(df, wave_cols, id_cols,
                              flat_range_threshold = 0.03,
                              low_thresh = 0.05, high_thresh = 0.95) {
  wave_mat <- as.matrix(df[, wave_cols])
  row_median <- apply(wave_mat, 1, median, na.rm = TRUE)
  row_range <- apply(wave_mat, 1, function(x) diff(range(x, na.rm = TRUE)))

  flagged <- row_range < flat_range_threshold &
    (row_median < low_thresh | row_median > high_thresh)

  df[flagged, id_cols, drop = FALSE] %>%
    mutate(
      spectrum_median = row_median[flagged],
      spectrum_range = row_range[flagged],
      check = "flat_extreme_spectrum"
    )
}

# ---- 2. Physically-impossible reflectance values ---------------------------

#' Flag any spectrum containing a wavelength reading outside a physically
#' possible reflectance range (small tolerance for instrument noise around
#' the 0/1 bounds, not a strict [0,1]).
qc_spectral_out_of_range <- function(df, wave_cols, id_cols,
                                      lower = -0.02, upper = 1.02) {
  wave_mat <- as.matrix(df[, wave_cols])
  n_bad <- rowSums(wave_mat < lower | wave_mat > upper, na.rm = TRUE)
  worst <- apply(wave_mat, 1, function(x) {
    x_bad <- x[x < lower | x > upper]
    if (length(x_bad) == 0) NA_real_ else x_bad[which.max(abs(x_bad - 0.5))]
  })

  flagged <- n_bad > 0
  df[flagged, id_cols, drop = FALSE] %>%
    mutate(
      n_bands_out_of_range = n_bad[flagged],
      worst_value = worst[flagged],
      check = "out_of_range_reflectance"
    )
}

# ---- 3. Within-group shape outlier (spectral angle from group median) -----

#' Flag spectra that differ strongly in shape from their own group's
#' median spectrum, using spectral angle (the angle in wavelength-space
#' between a spectrum and the group median spectrum -- insensitive to
#' overall brightness offsets, sensitive to shape). Flags a robust
#' (MAD-based) outlier in angle within each group.
#'
#' Designed for thin groups (few replicates per sampling unit, e.g. 6) --
#' deliberately not a PCA-distance method, which needs more replicates per
#' group than this kind of study typically has.
#'
#' @param group_var column defining "same expected spectrum" (e.g. a
#'   combined sample x treatment column)
qc_spectral_group_outliers <- function(df, group_var, wave_cols, id_cols,
                                        min_n = 4, z_threshold = 3) {
  wave_mat <- as.matrix(df[, wave_cols])

  spectral_angle <- function(a, b) {
    denom <- sqrt(sum(a^2)) * sqrt(sum(b^2))
    if (denom == 0) return(NA_real_)
    acos(pmin(1, pmax(-1, sum(a * b) / denom)))
  }

  df <- df %>% mutate(.row_id = row_number())

  results <- df %>%
    group_by(.data[[group_var]]) %>%
    filter(n() >= min_n) %>%
    group_modify(~ {
      idx <- .x$.row_id
      sub_mat <- wave_mat[idx, , drop = FALSE]
      group_median_spec <- apply(sub_mat, 2, median, na.rm = TRUE)
      angles <- apply(sub_mat, 1, spectral_angle, b = group_median_spec)

      med_angle <- median(angles, na.rm = TRUE)
      mad_angle <- mad(angles, na.rm = TRUE)
      robust_z <- if (mad_angle > 0) (angles - med_angle) / (1.4826 * mad_angle) else rep(0, length(angles))

      .x %>% mutate(spectral_angle_rad = angles, angle_robust_z = robust_z)
    }) %>%
    ungroup() %>%
    filter(abs(angle_robust_z) > z_threshold)

  if (nrow(results) == 0) {
    return(df[0, c(id_cols, group_var)] %>%
             mutate(spectral_angle_rad = numeric(0), angle_robust_z = numeric(0), check = character(0)))
  }

  results %>%
    select(all_of(id_cols), all_of(group_var), spectral_angle_rad, angle_robust_z) %>%
    mutate(check = "spectral_group_outlier") %>%
    arrange(desc(abs(angle_robust_z)))
}
