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
#' The reference is always a REAL measured scan from the group, never a
#' synthetic per-wavelength median (Henry, 2026-09-23: a per-wavelength
#' median is a spectrum nobody measured, which makes a flag against it hard
#' to interpret and hard to go back and look at). Two ways to pick it:
#'
#'   "medoid" - the scan with the smallest summed angle to all the others.
#'      Best in controlled lab settings (dry powder, press apparatus), where
#'      the scans in a group should be near-identical and the question is
#'      "which one doesn't belong".
#'   "albedo" - the scan whose albedo (mean reflectance across wave_cols) is
#'      closest to the group's median albedo. Expected to matter more for
#'      fresh-leaf and field spectroscopy, where overall brightness varies
#'      for real physical reasons.
#'
#' Set this per project in the config; Mill Test and other dry-powder
#' deposits use "medoid". The chosen method and the id of the scan actually
#' used are both returned, so they land in the QC log as provenance.
#'
#' The reference scan's own angle is 0 by construction, so it is excluded
#' from the median/MAD that define the z-scores (otherwise it drags the
#' centre down and inflates everyone else's z) and it is never itself
#' flagged.
#'
#' Limitation worth knowing: with 6 replicates the MAD is computed on 5
#' angles, which is a coarse scale estimate, so individual z values can be
#' large without meaning much. That is the price of working within such a
#' thin group, and it is exactly why qc_spectral_pca_flags() pools its
#' distance distribution globally instead. Run both and compare.
#'
#' Two ways to decide what counts as an outlier, set by `method`:
#'
#'   "absolute" (default) - flag when the angle to the reference exceeds
#'      `angle_max_deg`. Recommended. The angle is a physically meaningful
#'      quantity on a fixed scale, so a threshold in degrees means the same
#'      thing across projects and across reruns, and it is something a data
#'      user can actually interpret.
#'   "robust_z" - flag when the angle is a robust-z outlier among the other
#'      angles in its own group. Adapts to how tightly a given group happens
#'      to agree, but at typical replicate counts that is a liability rather
#'      than a feature: with 6 scans the MAD is computed on 5 angles, so the
#'      scale estimate is mostly noise, and a group whose scans agree
#'      unusually tightly will flag a perfectly good sibling for being
#'      slightly less tight. Kept for comparison, not recommended as the
#'      primary test.
#'
#' Either way the test is one-sided: only scans FURTHER from the reference
#' than expected are flagged. A small angle means good agreement.
#'
#' @param group_var column defining "same expected spectrum" (e.g. a
#'   combined sample x treatment column)
#' @param reference "medoid" or "albedo" -- see above
#' @param method "absolute" or "robust_z" -- see above
#' @param angle_max_deg absolute angle threshold in DEGREES, used when
#'   method = "absolute". Set this per project from the observed
#'   within-group distribution, not from a default.
qc_spectral_group_outliers <- function(df, group_var, wave_cols, id_cols,
                                        min_n = 4, z_threshold = 3,
                                        reference = c("medoid", "albedo"),
                                        method = c("absolute", "robust_z"),
                                        angle_max_deg = 1) {
  reference <- match.arg(reference)
  method <- match.arg(method)
  wave_mat <- as.matrix(df[, wave_cols])

  spectral_angle <- function(a, b) {
    denom <- sqrt(sum(a^2)) * sqrt(sum(b^2))
    if (denom == 0) return(NA_real_)
    acos(pmin(1, pmax(-1, sum(a * b) / denom)))
  }

  # Pairwise angles within a group, used to pick the medoid.
  angle_matrix <- function(m) {
    n <- nrow(m)
    out <- matrix(0, n, n)
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (i < j) {
          a <- spectral_angle(m[i, ], m[j, ])
          out[i, j] <- a; out[j, i] <- a
        }
      }
    }
    out
  }

  df <- df %>% mutate(.row_id = row_number())

  results <- df %>%
    group_by(.data[[group_var]]) %>%
    filter(n() >= min_n) %>%
    group_modify(~ {
      idx <- .x$.row_id
      sub_mat <- wave_mat[idx, , drop = FALSE]

      ref_i <- if (reference == "medoid") {
        which.min(rowSums(angle_matrix(sub_mat)))
      } else {
        albedo <- rowMeans(sub_mat, na.rm = TRUE)
        # With an even n the median falls between two scans; taking the
        # closest real scan picks one of them deterministically (ties -> first).
        which.min(abs(albedo - median(albedo, na.rm = TRUE)))
      }

      angles <- apply(sub_mat, 1, spectral_angle, b = sub_mat[ref_i, ])

      # Exclude the reference's own zero angle from the robust centre/scale.
      others <- angles[-ref_i]
      med_angle <- median(others, na.rm = TRUE)
      mad_angle <- mad(others, na.rm = TRUE)
      robust_z <- if (mad_angle > 0) (angles - med_angle) / (1.4826 * mad_angle) else rep(0, length(angles))
      robust_z[ref_i] <- 0

      .x %>% mutate(
        spectral_angle_deg = angles * 180 / pi,
        angle_robust_z = robust_z,
        is_reference = seq_len(nrow(.x)) == ref_i,
        reference_method = reference,
        reference_spectrum = .x[[id_cols[1]]][ref_i]
      )
    }) %>%
    ungroup() %>%
    # One-sided either way: a SMALL angle means the scan agrees well with the
    # reference, which is never a defect.
    filter(!is_reference,
           if (method == "absolute") spectral_angle_deg > angle_max_deg
           else angle_robust_z > z_threshold)

  if (nrow(results) == 0) {
    return(df[0, c(id_cols, group_var)] %>%
             mutate(spectral_angle_deg = numeric(0), angle_robust_z = numeric(0),
                    reference_method = character(0), reference_spectrum = character(0),
                    check = character(0)))
  }

  results %>%
    select(all_of(id_cols), all_of(group_var), spectral_angle_deg, angle_robust_z,
           reference_method, reference_spectrum) %>%
    mutate(check = "spectral_group_outlier",
           flag_method = method,
           threshold = if (method == "absolute") angle_max_deg else z_threshold) %>%
    arrange(desc(spectral_angle_deg))
}

# ---- 4. Single-band / band-range threshold --------------------------------

#' Flag spectra whose reflectance at a given band (or averaged over a band
#' range) falls outside a plausible range.
#'
#' Motivating case: the SHIFT dried-ground reflectance DAAC product removes
#' every spectrum with R350 >= 0.9 as "abnormal reflectance". Dry leaf
#' powder sits far below that at 350 nm, so a value that high means the
#' white reference panel or an empty puck got measured instead of the
#' sample. More direct than qc_spectral_flat() for that specific failure --
#' keep both, they catch different things.
#'
#' @param band a single wavelength, or c(lower, upper) to average over an
#'   inclusive range. Averaging 350-400 is worth considering because 350 nm
#'   on its own is the noisiest band on an ASD.
#' @param upper,lower bounds on the band value; supply either or both
#' @param wave_names numeric wavelengths matching wave_cols, in the same
#'   order. Defaults to parsing digits out of wave_cols (handles "X350").
qc_spectral_band_threshold <- function(df, wave_cols, id_cols, band,
                                        upper = NULL, lower = NULL,
                                        wave_names = NULL) {
  if (is.null(upper) && is.null(lower)) {
    stop("qc_spectral_band_threshold() needs at least one of upper/lower")
  }
  if (is.null(wave_names)) {
    wave_names <- as.numeric(gsub("[^0-9.]", "", wave_cols))
  }

  sel <- if (length(band) == 1) {
    which.min(abs(wave_names - band))
  } else {
    which(wave_names >= band[1] & wave_names <= band[2])
  }
  if (length(sel) == 0) stop("qc_spectral_band_threshold(): no wavelengths matched `band`")

  band_label <- if (length(band) == 1) {
    paste0("R", round(wave_names[sel]))
  } else {
    paste0("R", band[1], "-", band[2], "_mean")
  }

  vals <- if (length(sel) == 1) {
    as.numeric(df[[wave_cols[sel]]])
  } else {
    rowMeans(as.matrix(df[, wave_cols[sel], drop = FALSE]), na.rm = TRUE)
  }

  bad <- !is.na(vals) &
    ((!is.null(upper) & vals > upper) | (!is.null(lower) & vals < lower))

  df[bad, id_cols, drop = FALSE] %>%
    mutate(
      band = band_label,
      value = vals[bad],
      lower_bound = if (is.null(lower)) NA_real_ else lower,
      upper_bound = if (is.null(upper)) NA_real_ else upper,
      check = "spectral_band_threshold"
    )
}

# ---- 5. PCA review flags (two levels, one global fit) ----------------------

#' Flag spectra for MANUAL REVIEW using distances in a PCA space fit once
#' across the whole dataset (Henry, 2026-09-23: PCA flags data for review,
#' it never removes it automatically).
#'
#' Two levels:
#'   "within_sample"  - how far each replicate scan sits from its own
#'      sample's centroid. Catches one bad scan among the replicates.
#'   "between_sample" - how far each sample's centroid sits from the centroid
#'      of its group (species / treatment). Catches a whole sample that
#'      doesn't belong.
#'
#' The PCA is fit ONCE on every spectrum and both levels are measured in
#' that shared space. This is what makes the check work on thin replicate
#' sets: a per-group PCA on 6 scans is badly underpowered, which is why
#' qc_spectral_group_outliers() uses spectral angle instead. The two checks
#' are complementary -- run both.
#'
#' Distances are computed on the top `n_pc` components after scaling each
#' component by its global MAD (a diagonal robust Mahalanobis), so a
#' component with little spread doesn't get swamped by PC1. Because a single
#' group is far too small to robust-z within, the z-scores are taken across
#' the POOLED distance distribution -- every replicate for level 1, every
#' sample for level 2.
#'
#' @param sample_var column identifying one sample (the replicate set)
#' @param group_var column identifying the set a sample should resemble
#'   (species, treatment, species x treatment)
#' @param n_pc how many components to keep
qc_spectral_pca_flags <- function(df, wave_cols, id_cols, sample_var, group_var,
                                   n_pc = 5, within_z = 3, between_z = 3,
                                   min_n_within = 3, min_n_between = 3) {
  wave_mat <- as.matrix(df[, wave_cols])

  pca <- prcomp(wave_mat, center = TRUE, scale. = FALSE)
  n_pc <- min(n_pc, ncol(pca$x))
  scores <- pca$x[, seq_len(n_pc), drop = FALSE]

  pc_mad <- apply(scores, 2, mad, na.rm = TRUE)
  pc_mad[pc_mad == 0 | is.na(pc_mad)] <- 1  # degenerate component: leave unscaled
  scaled <- sweep(scores, 2, pc_mad, "/")

  var_explained <- sum(pca$sdev[seq_len(n_pc)]^2) / sum(pca$sdev^2)

  work <- df %>%
    select(all_of(unique(c(id_cols, sample_var, group_var)))) %>%
    mutate(.i = row_number())

  robust_z <- function(x) {
    m <- median(x, na.rm = TRUE); s <- mad(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    (x - m) / (1.4826 * s)
  }

  # --- level 1: replicate vs. its own sample centroid ---
  within <- work %>%
    group_by(.data[[sample_var]]) %>%
    filter(n() >= min_n_within) %>%
    group_modify(~ {
      s <- scaled[.x$.i, , drop = FALSE]
      centroid <- apply(s, 2, median, na.rm = TRUE)
      .x$pca_distance <- sqrt(rowSums(sweep(s, 2, centroid, "-")^2))
      .x
    }) %>%
    ungroup() %>%
    mutate(distance_robust_z = robust_z(pca_distance), level = "within_sample") %>%
    filter(distance_robust_z > within_z)

  # --- level 2: sample centroid vs. its group's centroid ---
  sample_centroids <- work %>%
    group_by(.data[[sample_var]], .data[[group_var]]) %>%
    summarise(across(everything(), ~ NA), .groups = "drop") %>%
    select(all_of(c(sample_var, group_var)))

  centroid_mat <- work %>%
    group_by(.data[[sample_var]]) %>%
    group_map(~ apply(scaled[.x$.i, , drop = FALSE], 2, median, na.rm = TRUE)) %>%
    do.call(rbind, .)
  centroid_keys <- work %>%
    group_by(.data[[sample_var]]) %>%
    group_keys()

  centroid_df <- bind_cols(centroid_keys, as.data.frame(centroid_mat)) %>%
    left_join(sample_centroids, by = sample_var)
  pc_names <- setdiff(names(centroid_df), c(sample_var, group_var))

  between <- centroid_df %>%
    group_by(.data[[group_var]]) %>%
    filter(n() >= min_n_between) %>%
    group_modify(~ {
      s <- as.matrix(.x[, pc_names, drop = FALSE])
      centroid <- apply(s, 2, median, na.rm = TRUE)
      .x$pca_distance <- sqrt(rowSums(sweep(s, 2, centroid, "-")^2))
      .x
    }) %>%
    ungroup() %>%
    mutate(distance_robust_z = robust_z(pca_distance), level = "between_sample") %>%
    filter(distance_robust_z > between_z) %>%
    select(all_of(c(sample_var, group_var)), pca_distance, distance_robust_z, level)

  out <- bind_rows(
    within %>% select(all_of(id_cols), any_of(c(sample_var, group_var)),
                      pca_distance, distance_robust_z, level),
    between
  )

  if (nrow(out) == 0) {
    return(tibble::tibble(pca_distance = numeric(0), distance_robust_z = numeric(0),
                          level = character(0), check = character(0)))
  }

  out %>%
    mutate(check = "spectral_pca_review",
           n_pc = n_pc,
           pc_var_explained = round(var_explained, 4)) %>%
    arrange(desc(distance_robust_z))
}

# ---- 6. Splice / detector-boundary jump --------------------------------

#' Flag unusually large steps at an instrument's detector splice points,
#' measured BEFORE jump correction.
#'
#' Jump correction hides fiber, probe and detector-temperature problems by
#' construction -- once the segments have been made continuous, the evidence
#' is gone. So this has to run on genuinely uncorrected spectra. Worth
#' knowing: a file named "raw" is not necessarily uncorrected. A corrected
#' spectrum typically has R(splice) == R(splice+1) exactly, and those may be
#' the only exactly-equal adjacent pairs in the whole spectrum -- a cheap
#' way to tell (see qc_detect_jump_corrected()).
#'
#' Matters most for field spectra, where detector temperature shifts the
#' jumps and the manufacturer's temperature-correction curve usually isn't
#' available. Lower priority for dry-powder lab data.
#'
#' @param splices wavelengths of the detector boundaries. Supply from the
#'   per-project config -- see qc_splice_defaults().
#' @param jump_threshold flag when |R(splice+1) - R(splice)| exceeds this
qc_spectral_splice_jump <- function(df, wave_cols, id_cols, splices,
                                     jump_threshold = 0.02, wave_names = NULL) {
  if (is.null(wave_names)) wave_names <- as.numeric(gsub("[^0-9.]", "", wave_cols))
  wave_mat <- as.matrix(df[, wave_cols])

  purrr::map_dfr(splices, function(sp) {
    lo_i <- which.min(abs(wave_names - sp))
    hi_i <- lo_i + 1
    if (hi_i > ncol(wave_mat)) return(NULL)

    jump <- wave_mat[, hi_i] - wave_mat[, lo_i]
    bad <- !is.na(jump) & abs(jump) > jump_threshold
    if (!any(bad)) return(NULL)

    df[bad, id_cols, drop = FALSE] %>%
      mutate(splice_nm = wave_names[lo_i],
             jump = jump[bad],
             jump_threshold = jump_threshold,
             check = "spectral_splice_jump")
  })
}

#' Cheap test for whether spectra have already been jump-corrected: a
#' boundary-matching correction leaves R(splice) == R(splice+1) exactly,
#' which is otherwise vanishingly rare in real data. Returns one row per
#' splice with the share of spectra showing an exact match.
qc_detect_jump_corrected <- function(df, wave_cols, splices, wave_names = NULL) {
  if (is.null(wave_names)) wave_names <- as.numeric(gsub("[^0-9.]", "", wave_cols))
  wave_mat <- as.matrix(df[, wave_cols])

  purrr::map_dfr(splices, function(sp) {
    lo_i <- which.min(abs(wave_names - sp))
    hi_i <- lo_i + 1
    if (hi_i > ncol(wave_mat)) return(NULL)
    exact <- wave_mat[, hi_i] == wave_mat[, lo_i]
    tibble::tibble(
      splice_nm = wave_names[lo_i],
      n_spectra = nrow(wave_mat),
      prop_exactly_equal = mean(exact, na.rm = TRUE),
      likely_corrected = mean(exact, na.rm = TRUE) > 0.9
    )
  })
}

#' Known detector splice wavelengths by instrument.
#'
#' Both sets below are taken from the lab's own jump-correction code
#' (github.com/krkovach/SpectralPredictR, main branch, `apply_jump_correction`,
#' as vendored into Mill Test's SpectralPredictR.R) -- ASD c(1000, 1800) and
#' .sig/SVC c(990, 1900). They are not independent readings of the
#' manufacturer documentation; if a project needs them verified against the
#' instrument manual, do that rather than trusting this function.
#'
#' Accepts instrument names case-insensitively and treats "svc" and "sig" as
#' the same instrument (the SVC writes .sig files). This alias is the point:
#' upstream `apply_jump_correction()` lowercases its argument and then tests
#' only `== "sig"`, so `instrument = "SVC"` matches nothing, leaves `splices`
#' NULL, and the correction then misbehaves. Harmless for ASD-only work,
#' broken for SVC.
#'
#' PSR (Spectral Evolution) is deliberately absent rather than guessed, per
#' Henry's instruction (2026-09-23) to take PSR detector boundaries from the
#' instrument documentation. Add them here once that's in hand.
qc_splice_defaults <- function(instrument) {
  key <- tolower(instrument)
  if (key %in% c("svc", "sig")) key <- "svc"
  switch(key,
    asd = c(1000, 1800),
    svc = c(990, 1900),
    stop(sprintf(
      "No splice defaults for instrument '%s'. Supply `splices` explicitly in the project config. (PSR boundaries are intentionally not hard-coded -- take them from the instrument documentation.)",
      instrument))
  )
}
