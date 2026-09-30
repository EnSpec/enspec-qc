##################################################
#
# spectral_transforms.R
#
# Per-spectrum transforms applied before the shape-based spectral QC checks
# (qc_spectral_group_outliers(), qc_spectral_pca_flags()). The physical
# checks (flat, out of range, band thresholds, splice jumps) must still run
# on reflectance: a transform that removes offset and scale also removes the
# evidence those checks look for.
#
# snv(): standard normal variate. Barnes, Dhanoa & Lister (1989) Appl.
# Spectrosc. 43:772-777, eq. 1: center each spectrum on its own mean and
# divide by its own standard deviation (n - 1). Removes additive offset and
# multiplicative scatter, i.e. most particle-size and packing differences
# between wells, so the shape checks compare spectral shape (Henry,
# 2026-09-30). In the Mill Test grind experiment SNV cut the grind share of
# spectral variance from 23% to 5% (EnSpec/dryspec-preprocessing,
# 04_mill24_mechanism.R).
#
# Moved here 2026-09-30 from ACRES (workflows/5_daac_release_prep/code/
# spectral_transforms.R), which copied snv() from EnSpec/dryspec-preprocessing
# R/transforms.R, where it is tested against prospectr::standardNormalVariate.
# Pure functions, no paths.
#
##################################################

snv <- function(X) {
  X <- as.matrix(X)
  mu <- rowMeans(X)
  s  <- sqrt(rowSums((X - mu)^2) / (ncol(X) - 1))
  (X - mu) / s
}

#' SNV over the wavelength columns within `range` (inclusive); returns a
#' data frame with `keep_cols` plus the transformed columns only. The SNV
#' mean and SD are computed over `range`, so a different range gives a
#' different transform; record the range in the project config.
snv_over_range <- function(df, wave_cols, range, keep_cols) {
  wl <- as.numeric(gsub("[^0-9.]", "", wave_cols))
  sel <- wave_cols[wl >= range[1] & wl <= range[2]]
  if (length(sel) < 2) {
    stop(sprintf("snv_over_range(): %d wavelength columns fall in %g-%g nm; check wave_cols and range",
                 length(sel), range[1], range[2]))
  }
  out <- as.data.frame(snv(df[, sel, drop = FALSE]))
  names(out) <- sel
  cbind(df[, keep_cols, drop = FALSE], out)
}
