##################################################
#
# qc_report.R
#
# Provenance/logging helpers shared across QC checks. Every check should
# write its flags through these, so a deposit always ships: cleaned data,
# flag columns carried into the deposit, and a short QC report -- nothing
# dropped silently (per the 2026-09-20 QC standard decision).
#
##################################################

library(dplyr)

#' Write one check's flag table to disk with a consistent naming
#' convention, and return it invisibly (so this can sit inline in a pipe).
qc_write_log <- function(flag_df, out_dir, check_name) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, paste0(check_name, "_flags.csv"))
  write.csv(flag_df, path, row.names = FALSE)
  message(sprintf("  %s: %d rows flagged -> %s", check_name, nrow(flag_df), path))
  invisible(flag_df)
}

#' Summarize a named list of per-check flag tables into one counts table.
#'
#' Deliberately stops here rather than also building a per-id combined
#' `quality_flag` column: checks in this library are allowed to flag at
#' different grains (a per-replicate spectrum table vs. a per-sample trait
#' table vs. a per-sample x treatment table), and only the calling project
#' knows how its own id columns relate to each other -- e.g. broadcasting a
#' sample-level flag down to every spectrum row for that sample is a join
#' the project has to define, not something this library can infer. Do that
#' join in the project's own script, then write the combined per-row
#' quality_flag column with qc_write_log() alongside everything else.
#'
#' @param flag_list named list of flag data frames (as returned by the
#'   qc_* functions in this package)
qc_build_report <- function(flag_list, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  check_counts <- bind_rows(lapply(names(flag_list), function(nm) {
    tibble::tibble(check_name = nm, n_flagged = nrow(flag_list[[nm]]))
  }))

  write.csv(check_counts, file.path(out_dir, "qc_check_summary.csv"), row.names = FALSE)

  message("\n== QC summary ==")
  print(check_counts)

  check_counts
}

#' Apply removals (NA-out) for checks tagged as "remove" policy in a
#' project's config, logging exactly what was changed. Flag-only checks
#' should never be passed here -- keep this call explicit and separate
#' from the flagging step so removal is always an opt-in, visible action.
#'
#' @param df data frame to modify
#' @param removal_flags long-format flag table (id_col, trait_flagged) for
#'   rows/traits to NA out
#' @param id_col column identifying rows in df (must match removal_flags)
#' @param trait_col column in removal_flags naming which trait to NA
#'   (defaults to "trait_flagged", the convention used by every qc_* check)
qc_apply_removals <- function(df, removal_flags, id_col, out_dir,
                              trait_col = "trait_flagged") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  removal_log <- list()

  for (tr in unique(removal_flags[[trait_col]])) {
    if (!tr %in% names(df)) next
    ids <- removal_flags[[id_col]][removal_flags[[trait_col]] == tr]
    affected <- df[[id_col]] %in% ids
    if (!any(affected)) next

    removal_log[[tr]] <- df[affected, id_col, drop = FALSE] %>%
      mutate(trait = tr, removed_value = df[[tr]][affected], removal_date = Sys.Date())

    df[[tr]][affected] <- NA
  }

  removal_log_df <- bind_rows(removal_log)
  write.csv(removal_log_df, file.path(out_dir, "removal_log.csv"), row.names = FALSE)
  message(sprintf("Removed (set NA) %d values across %d traits -> %s",
                  nrow(removal_log_df), length(removal_log), file.path(out_dir, "removal_log.csv")))

  list(data = df, removal_log = removal_log_df)
}
