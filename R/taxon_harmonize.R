##################################################
#
# taxon_harmonize.R
#
# Shared taxon-name harmonization. Generalized from the vetted GCFR pipeline
# (Data_Workflow_2_Taxonomic_Cleaning/Code/Workflow2_Taxonomic_Clean.R, Henry
# Frye), which is the most recently reviewed version of a workflow Henry has
# now written several times over.
#
# Runs BEFORE the QC checks in trait_qc.R / spectral_qc.R, or beside them: the
# QC grouping variable is usually taxon, so harmonization has to be settled
# first (Henry, 2026-09-23).
#
# Same split as the rest of this library. General here: the normalization
# rules, the match cascade, the join safety, the recovery classification.
# Per-project, in a config: the regional checklist, the botanist-verified
# override table, the placeholder patterns, the fuzzy threshold.
#
# Two hard-won behaviours from the GCFR work are built in as guards rather
# than left to the caller, because both failed silently there:
#
#   1. The WFO join MUST stay 1:1. GCFR's original join was on a natural key
#      that wasn't unique, which quietly turned 9,547 rows into 27,337 (and
#      2,509 into 7,431). taxon_match_wfo() joins on a synthetic row index and
#      asserts the row count, so that cannot recur.
#   2. Infraspecific epithets get dropped. A reviewer caught subsp./var.
#      disappearing from some records mid-pipeline; it took a dedicated
#      diagnostic to find where. taxon_check_infraspecific_loss() turns that
#      into a standing check.
#
# Nothing here drops a name. An unmatched name is classified, never discarded
# -- per the 2026-09-20 QC standard, nothing is removed silently.
#
##################################################

library(dplyr)
library(rlang)

# ---- 1. Name string normalization -----------------------------------------

#' Normalize taxon name strings: infraspecific abbreviations, spacing, and
#' plot-based placeholder names.
#'
#' These are the fixes that make otherwise-identical names fail to match each
#' other or the backbone: "var" without its period, "cf.crassa" with no space,
#' doubled spaces from an earlier edit.
#'
#' @param placeholder_patterns regexes for field placeholders that are not
#'   taxon names at all (e.g. "Plot70sp1"). Project-specific, since every
#'   field crew invents its own. Matches are replaced with `placeholder_label`
#'   so they are explicitly unknown rather than silently fuzzy-matched to some
#'   unrelated species.
taxon_standardize_names <- function(x,
                                    placeholder_patterns = "^Plot\\d+sp\\d*$",
                                    placeholder_label = "unknown species",
                                    infra_ranks = c("subsp", "var", "cf")) {
  ranks <- paste(infra_ranks, collapse = "|")

  # add a period after subsp/var/cf when missing, without doubling an existing one
  x <- gsub(sprintf("\\b(%s)\\b(?!\\.)", ranks), "\\1.", x, perl = TRUE)
  # exactly one space after the abbreviation ("cf.crassa" -> "cf. crassa")
  x <- gsub(sprintf("\\b(%s)\\.\\s*", ranks), "\\1. ", x, perl = TRUE)
  # collapse whitespace introduced above or already present
  x <- trimws(gsub("\\s+", " ", x))

  for (pat in placeholder_patterns) {
    x <- ifelse(grepl(pat, x, ignore.case = TRUE), placeholder_label, x)
  }
  x
}

#' Strip indeterminate suffixes ("sp", "sp.", "sp 1", "species") so a
#' genus-only record resolves to its genus.
#'
#' Critically, this must NOT strip subsp./var./cf. -- those are real
#' infraspecific information, and "subsp" contains "sp". The negative check
#' below is the whole reason this is a named function rather than an inline
#' gsub: the GCFR releve data went unresolved through two separate pipeline
#' passes because this step had been applied to some inputs and not others.
taxon_strip_indet <- function(x, infra_ranks = c("subsp", "var", "cf")) {
  ranks <- paste(infra_ranks, collapse = "|")
  indet <- "\\bsp\\.?\\b|\\bsp\\s*\\d+\\b|\\bspecies\\b"

  ifelse(
    grepl(indet, x, ignore.case = TRUE) &
      !grepl(sprintf("\\b(%s)\\b", ranks), x, ignore.case = TRUE),
    gsub("\\s+(sp\\.?\\s*\\d*|species).*", "", x, ignore.case = TRUE),
    x
  )
}

#' Parse a name into genus / species epithet / infraspecific rank + epithet.
taxon_parse_name <- function(x, infra_ranks = c("subsp", "var", "cf")) {
  ranks <- paste(infra_ranks, collapse = "|")
  x <- stringr::str_squish(as.character(x))
  tibble::tibble(
    raw             = x,
    genus           = stringr::word(x, 1),
    species_epithet = stringr::word(x, 2),
    infra_rank      = stringr::str_extract(x, sprintf("\\b(%s)\\.?\\b", ranks)),
    infra_epithet   = stringr::str_extract(x, sprintf("(?<=(%s)\\.?\\s)\\S+", ranks))
  )
}

# ---- 2. Botanist-verified manual overrides --------------------------------

#' Apply manual name corrections, and report which ones actually fired.
#'
#' Some names no automated matcher will get right, and a botanist's call is
#' the authority -- in GCFR, WFO resolved `Erica demissa` to the wrong
#' accepted name, and a narrow endemic had been recorded well outside its
#' range. Those corrections belong in a reviewable table, not scattered
#' through a script as `case_when` arms.
#'
#' Reporting which overrides fired matters as much as applying them: an
#' override that no longer matches anything is either a fixed upstream error
#' or a silently broken rule, and both are worth seeing.
#'
#' @param overrides a data frame with `from` and `to` columns, optionally
#'   `record_id` to scope a correction to one record (the same misspelling may
#'   be right elsewhere), and optionally `reason` for provenance.
#' @param id_col the column in `df` that `overrides$record_id` refers to.
#' @return `df` with `name_col` corrected, plus an `overrides_applied`
#'   attribute holding the per-rule hit counts.
taxon_apply_overrides <- function(df, name_col, overrides, id_col = NULL) {
  stopifnot(all(c("from", "to") %in% names(overrides)))
  scoped <- "record_id" %in% names(overrides) && !is.null(id_col)

  hits <- integer(nrow(overrides))
  out <- df[[name_col]]

  for (i in seq_len(nrow(overrides))) {
    match_i <- out == overrides$from[i]
    if (scoped && !is.na(overrides$record_id[i])) {
      match_i <- match_i & df[[id_col]] == overrides$record_id[i]
    }
    match_i[is.na(match_i)] <- FALSE
    hits[i] <- sum(match_i)
    out[match_i] <- overrides$to[i]
  }

  df[[name_col]] <- out

  applied <- overrides %>%
    mutate(n_records_changed = hits,
           fired = hits > 0)

  unused <- applied %>% filter(!fired)
  if (nrow(unused) > 0) {
    warning(sprintf("%d taxon override(s) matched nothing: %s. Either the upstream data was fixed (remove the rule) or the rule no longer matches (check it).",
                    nrow(unused), paste(unused$from, collapse = ", ")),
            call. = FALSE)
  }

  attr(df, "overrides_applied") <- applied
  df
}

# ---- 3. Comparison against a regional checklist ---------------------------

#' Report names absent from a regional checklist / taxonomic authority.
#'
#' Flag-only and informational. A regional flora (Goldblatt & Manning for the
#' Cape, say) is the right authority for its own region but is not a global
#' one, so absence means "look at this", never "wrong". Sarah's regional-flora
#' work is exactly this case.
taxon_check_against_checklist <- function(df, name_col, checklist,
                                          checklist_col = "Taxon",
                                          id_cols = NULL) {
  accepted <- unique(checklist[[checklist_col]])

  df %>%
    filter(!is.na(.data[[name_col]]), !.data[[name_col]] %in% accepted) %>%
    select(all_of(c(id_cols, name_col))) %>%
    mutate(check = "taxon_not_in_checklist",
           checklist_size = length(accepted)) %>%
    arrange(.data[[name_col]])
}

# ---- 4. WFO matching, with the 1:1 join guaranteed ------------------------

#' Resolve names against a WFO static backbone, keeping the join strictly 1:1.
#'
#' `WFO.match()` / `WFO.one()` return exactly one row per input row, so a
#' synthetic row index is a safe join key -- and a natural key is not, unless
#' you have proven it unique. GCFR joined on `NewUID`, which needed collector
#' and replicate to be unique, and the result was a many-to-many blow-up
#' (9,547 rows -> 27,337) that was not noticed until a later workflow. The
#' assertions below are the point of this function.
#'
#' @param wfo_data a loaded WFO backbone (see wfo_backbone_check() first --
#'   record which version produced these names)
#' @param fuzzy passed to WFO.match's Fuzzy.min
#' @param keep WFO columns to carry back
taxon_match_wfo <- function(df, name_col, wfo_data, fuzzy = TRUE,
                            keep = c("scientificName", "family",
                                     "scientificNameAuthorship")) {
  if (!requireNamespace("WorldFlora", quietly = TRUE)) {
    stop("taxon_match_wfo() requires the 'WorldFlora' package")
  }
  if (".wfo_row_id" %in% names(df)) {
    stop("df already has a .wfo_row_id column; rename it before calling")
  }

  df$.wfo_row_id <- seq_len(nrow(df))

  matched <- WorldFlora::WFO.match(spec.data = df, WFO.data = wfo_data,
                                   spec.name = name_col, Fuzzy.min = fuzzy)
  single <- WorldFlora::WFO.one(matched)

  taxon_join_one_to_one(df, single, keep = keep)
}

#' Join a matcher's per-row output back onto its input, refusing to proceed
#' unless the relationship really is 1:1.
#'
#' Split out from taxon_match_wfo() so this -- the part that actually prevents
#' the GCFR row blow-up -- can be tested without a 950 MB backbone or the
#' WorldFlora package, and reused for any matcher with the same one-row-per-
#' input contract (GBIF, TNRS, a regional flora lookup).
#'
#' Fails loudly rather than returning a silently multiplied table, because
#' that is precisely how the original bug survived into a later workflow.
#'
#' @param df the input, carrying `key`
#' @param matched the matcher's output, expected one row per input row
#' @param key the synthetic row index both sides share
taxon_join_one_to_one <- function(df, matched, keep, key = ".wfo_row_id") {
  if (!key %in% names(matched)) {
    stop(sprintf("matcher output has no '%s' column, so the join cannot be verified 1:1. Do not fall back to a natural key unless you have proven it unique.",
                 key))
  }
  if (nrow(matched) != nrow(df)) {
    stop(sprintf("matcher returned %d rows for %d inputs -- the 1:1 contract this join relies on does not hold; do not trust the result.",
                 nrow(matched), nrow(df)))
  }
  if (any(duplicated(matched[[key]]))) {
    stop(sprintf("matcher output has duplicate '%s' values -- joining on it would multiply rows.", key))
  }

  keep <- intersect(keep, names(matched))
  out <- left_join(df, matched[, c(key, keep), drop = FALSE], by = key)

  if (nrow(out) != nrow(df)) {
    stop(sprintf("the join changed the row count (%d -> %d); it was not 1:1.",
                 nrow(df), nrow(out)))
  }

  out[, setdiff(names(out), key), drop = FALSE]
}

# ---- 5. Recovery classification ------------------------------------------

#' Classify every record by how far its name got resolved.
#'
#' The reason this exists rather than just leaving NAs: an unresolved name and
#' a name resolved only to family and a non-taxonomic junk entry are three
#' different data-quality situations, and a downstream user needs to tell them
#' apart. A single NA collapses them.
#'
#' @param matched_col the resolved scientific name (NA where unmatched)
#' @param family_col optional family, which may be recoverable even when the
#'   name is not (a record given only as "Orchid", say)
#' @param non_taxonomic a predicate (logical vector) marking entries that were
#'   never taxon names -- stray numbers, codes, notes
taxon_classify_recovery <- function(df, matched_col, family_col = NULL,
                                    non_taxonomic = NULL) {
  has_name <- !is.na(df[[matched_col]])
  has_family <- if (is.null(family_col)) rep(FALSE, nrow(df)) else !is.na(df[[family_col]])
  junk <- if (is.null(non_taxonomic)) rep(FALSE, nrow(df)) else non_taxonomic

  df %>%
    mutate(recovery_status = case_when(
      has_name    ~ "wfo_matched",
      junk        ~ "flagged_non_taxonomic_entry",
      has_family  ~ "recovered_family_only",
      TRUE        ~ "unresolved_morphotype"
    ))
}

# ---- 6. Infraspecific-epithet loss check ---------------------------------

#' Flag records that had an infraspecific epithet before harmonization and
#' lost it after.
#'
#' A standing guard for the failure a GCFR reviewer found: subsp./var.
#' silently disappearing from some records and some columns but not others,
#' which needed a bespoke diagnostic to locate. Losing one is not always wrong
#' -- WFO may legitimately resolve an infraspecific name to the species -- but
#' it should be a visible, counted decision rather than a surprise at review.
#'
#' @param before,after vectors of the same length: names entering and leaving
#'   harmonization
taxon_check_infraspecific_loss <- function(df, before, after, id_cols = NULL,
                                           infra_ranks = c("subsp", "var", "cf")) {
  ranks <- sprintf("\\b(%s)\\b", paste(infra_ranks, collapse = "|"))

  df %>%
    mutate(
      .name_before = .data[[before]],
      .name_after  = .data[[after]],
      had_infra  = grepl(ranks, .name_before, ignore.case = TRUE),
      has_infra  = grepl(ranks, .name_after, ignore.case = TRUE) &
                     !is.na(.name_after)
    ) %>%
    filter(had_infra, !has_infra) %>%
    select(all_of(id_cols), name_before = .name_before, name_after = .name_after) %>%
    mutate(check = "infraspecific_epithet_lost",
           reason = "name carried subsp./var./cf. before harmonization and not after -- may be a legitimate resolution to species, but confirm it is intended") %>%
    distinct()
}
