##################################################
#
# try_reference.R
#
# Helper for pulling literature-range bounds out of a bulk TRY database
# export, for use as qc_hard_bounds() / flag-only literature-range bounds.
# Generalizes the pattern from the GCFR trait QC (Workflow3_TRY_data_ranges.R).
#
##################################################

#' Compute min/max ranges for a set of TRY trait names from a raw TRY
#' database text export (tab-separated, one row per trait record).
#'
#' @param try_path path to the raw TRY .txt export
#' @param trait_names character vector of TraitName values to pull (must
#'   match TRY's own naming exactly -- check with the commented-out full
#'   listing code below if unsure)
#' @param error_risk_cutoff TRY's ErrorRisk is a standardized outlier score;
#'   TRY's own docs suggest >3-4 is often treated as suspect. No universal
#'   hard cutoff -- adjust per trait if needed.
#' @param signed_traits trait names allowed to be negative (e.g. isotope
#'   signatures) -- all others are filtered to StdValue > 0
qc_try_bounds <- function(try_path, trait_names,
                          error_risk_cutoff = 4,
                          signed_traits = character(0)) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("qc_try_bounds() requires the 'data.table' package")
  }
  try_data <- data.table::fread(try_path, sep = "\t", quote = "")

  subset_data <- try_data[
    TraitName %in% trait_names &
      !is.na(StdValue) &
      (is.na(ErrorRisk) | ErrorRisk <= error_risk_cutoff) &
      (TraitName %in% signed_traits | StdValue > 0)
  ]

  subset_data[, .(min_val = min(StdValue), max_val = max(StdValue),
                  n = .N, unit = data.table::first(UnitName)),
             by = TraitName]
}

# To see the full list of trait names available in a given TRY export:
#   dt <- data.table::fread(try_path, sep = "\t", quote = "", select = "TraitName")
#   sort(unique(dt$TraitName))
