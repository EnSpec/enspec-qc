##################################################
#
# wfo_backbone.R
#
# Version/currency checks for the World Flora Online static taxonomic
# backbone, the reference most of this lab's taxon harmonization resolves
# against.
#
# Why this exists (Henry, 2026-09-24): WFO publishes the backbone as a dated
# static file, and it is very easy to keep resolving names against whichever
# copy happens to be on the drive -- for years. The file works, the matching
# succeeds, nothing errors, and the names you get back are simply out of date.
# This came up in review.
#
# Design: OFFLINE-FIRST, deliberately. The strongest check needs no network at
# all -- if a newer backbone is already sitting in the same directory as the
# one you are using, that is the mistake, stated plainly. A second check asks
# only "how old is the file I am using", which flags a lab-wide stale copy
# without needing to know what the newest release is.
#
# There is no reliable machine-readable source for "the current WFO version".
# worldfloraonline.org/downloadData is HTML with no version API, and as of
# 2026-09-24 its TLS chain fails to verify from at least one lab machine
# (`curl` returns "unable to verify the first certificate"), so a scraper
# would be both untestable and a runtime dependency on a host we cannot
# reliably reach. Zenodo carries WFO's descriptive record but not the periodic
# versioned data releases. So instead of scraping, the caller can record the
# newest version they have actually seen as `latest_known` in the project
# config, which makes the check exact, and the warning always names the
# download page so a human can confirm.
#
##################################################

library(dplyr)

WFO_DOWNLOAD_URL <- "https://www.worldfloraonline.org/downloadData"

#' Parse a WFO backbone version out of a file path.
#'
#' WFO names these files `classification_v.YYYY.MM.csv` (e.g.
#' `classification_v.2026.6.csv`), with the month not zero-padded in at least
#' some releases, so both `2026.6` and `2026.06` are accepted.
#'
#' @return a one-row tibble: version, year, month, released (the first of that
#'   month, as a Date), and age_months relative to `as_of`. Returns NA fields
#'   rather than erroring if the path carries no recognizable version, because
#'   a caller may legitimately be using a renamed or custom backbone -- that
#'   case is reported by wfo_backbone_check() instead.
wfo_backbone_version <- function(path, as_of = Sys.Date()) {
  v <- stringr::str_match(basename(path), "v\\.?(\\d{4})\\.(\\d{1,2})")

  if (is.na(v[1, 1])) {
    return(tibble::tibble(path = path, version = NA_character_,
                          year = NA_integer_, month = NA_integer_,
                          released = as.Date(NA), age_months = NA_real_))
  }

  year  <- as.integer(v[1, 2])
  month <- as.integer(v[1, 3])
  released <- as.Date(sprintf("%04d-%02d-01", year, month))

  tibble::tibble(
    path = path,
    version = sprintf("%d.%d", year, month),
    year = year, month = month, released = released,
    age_months = as.numeric(difftime(as_of, released, units = "days")) / 30.44
  )
}

#' Check that a WFO backbone is the newest one available, and record which one
#' was used.
#'
#' Recording is not a secondary concern here: a harmonized name is only
#' reproducible if you know which backbone produced it, so this returns the
#' version for the caller to write into its provenance whether or not any
#' check trips.
#'
#' @param path the backbone file actually being used
#' @param latest_known optionally, the newest version the project knows to
#'   exist, as "YYYY.M" or "YYYY.MM". Set it in the project config and bump it
#'   when someone checks the download page; that makes this check exact
#'   instead of merely age-based.
#' @param max_age_months warn if the backbone in use is older than this.
#'   Default 18 months -- long enough not to nag about a deliberate pin,
#'   short enough that a genuinely abandoned copy surfaces.
#' @param search_dir also look for newer backbones alongside this one. This is
#'   the check that catches the real mistake, so it is on by default.
#' @return a list with `version` (a one-row tibble) and `findings` (a tibble,
#'   zero rows if the backbone is current). Warns once per finding.
wfo_backbone_check <- function(path, latest_known = NULL,
                               max_age_months = 18, search_dir = TRUE,
                               as_of = Sys.Date()) {
  in_use <- wfo_backbone_version(path, as_of = as_of)
  findings <- list()

  if (!file.exists(path)) {
    findings[["missing"]] <- tibble::tibble(
      issue = "backbone_not_found", severity = "error",
      detail = sprintf("no file at %s", path))
  }

  if (is.na(in_use$version)) {
    findings[["unparsed"]] <- tibble::tibble(
      issue = "version_unrecognized", severity = "warning",
      detail = sprintf("cannot read a WFO version from '%s'. Expected a name like classification_v.2026.6.csv. Keep the version in the filename -- a harmonized name is only reproducible if the backbone that produced it is identifiable.",
                       basename(path)))
  }

  # The one that matters: a newer backbone already on disk next to this one.
  if (search_dir && !is.na(in_use$version)) {
    siblings <- list.files(dirname(path), pattern = "v\\.?\\d{4}\\.\\d{1,2}",
                           full.names = TRUE)
    siblings <- setdiff(siblings, path)
    if (length(siblings) > 0) {
      sib <- bind_rows(lapply(siblings, wfo_backbone_version, as_of = as_of)) %>%
        filter(!is.na(released), released > in_use$released)
      if (nrow(sib) > 0) {
        newest <- sib %>% arrange(desc(released)) %>% slice(1)
        findings[["newer_local"]] <- tibble::tibble(
          issue = "newer_backbone_on_disk", severity = "warning",
          detail = sprintf("using v.%s but v.%s is already in the same directory (%s). This is the common mistake: the old file still works, so nothing complains.",
                           in_use$version, newest$version, basename(newest$path)))
      }
    }
  }

  if (!is.null(latest_known) && !is.na(in_use$version)) {
    latest <- wfo_backbone_version(paste0("v.", latest_known), as_of = as_of)
    if (!is.na(latest$released) && latest$released > in_use$released) {
      findings[["behind_known"]] <- tibble::tibble(
        issue = "newer_backbone_published", severity = "warning",
        detail = sprintf("using v.%s but the project config records v.%s as available. Download: %s",
                         in_use$version, latest$version, WFO_DOWNLOAD_URL))
    }
  }

  if (!is.na(in_use$age_months) && in_use$age_months > max_age_months) {
    findings[["stale"]] <- tibble::tibble(
      issue = "backbone_stale", severity = "warning",
      detail = sprintf("v.%s is %.0f months old (threshold %.0f). Check for a newer release: %s",
                       in_use$version, in_use$age_months, max_age_months,
                       WFO_DOWNLOAD_URL))
  }

  findings <- if (length(findings) == 0) {
    tibble::tibble(issue = character(0), severity = character(0),
                   detail = character(0))
  } else {
    bind_rows(findings)
  }

  if (nrow(findings) > 0) {
    for (i in seq_len(nrow(findings))) {
      warning(sprintf("WFO backbone: %s", findings$detail[i]), call. = FALSE)
    }
  } else {
    message(sprintf("WFO backbone v.%s (%s) -- current as far as this check can tell.",
                    in_use$version, basename(path)))
  }

  list(version = in_use, findings = findings)
}
