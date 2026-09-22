# enspec-qc

Shared, project-agnostic trait and spectral QC functions for EnSpec lab
data releases (field/lab trait data + proximal/lab spectroscopy).

**Status:** local only, not yet pushed to a remote. First real
implementation was built alongside the Mill Test 2024 EcoSIS deposit
(2026-09); adopted from there into other projects as they need it.

## Design

One shared, project-agnostic set of checks, with everything
project-specific supplied by the caller: the grouping variable that
defines "expected variability" (a taxon, a treatment, a taxon x treatment
combination, ...), trait bounds, and thresholds. The comparison logic here
is general; what counts as a sampling unit or a plausible value is not.

- **Trait QC** (`R/trait_qc.R`):
  1. `qc_hard_bounds()` -- flag values outside a physically-possible or
     literature-supported range (mis-entered data: decimal shifts, sign
     errors, unit mismatches).
  2. `qc_bad_single_measure()` / `qc_bad_ratio_component()` -- flag a
     replicate measurement (or a ratio trait's raw component) that's off
     by a large factor from the leave-one-out median of its own sampling
     unit. Catches balance/instrument entry errors, not real biological
     variation, so it operates on tight, same-individual/same-session
     groups (e.g. replicate weights for one sample).
  3. `qc_magnitude_shift()` -- the same idea at a coarser grain (a group's
     median, e.g. species or species x treatment), for catching
     decimal-shift errors that wouldn't show up in a small replicate set.
  4. `qc_distribution_outliers()` -- robust (MAD-based) z-score flag
     within a group. Flag-only: distribution outliers can be real biology.
  5. `qc_covariation_outliers()` -- robust SMA regression between two
     traits expected to strongly covary (e.g. N vs. LMA); flags points far
     from the fit. Flag-only.
  6. `qc_try_bounds()` (`R/try_reference.R`) -- pulls min/max literature
     ranges for a set of traits out of a bulk TRY database export, for use
     as `qc_hard_bounds()` input.

- **Spectral QC** (`R/spectral_qc.R`) -- for contact-probe / press-apparatus
  reflectance spectra (not field/canopy spectroscopy). No prior internal
  precedent; built fresh for Mill Test:
  1. `qc_spectral_flat()` -- flags spectra that are both flat (low
     variance across wavelength) and sitting at an extreme reflectance
     level (near 0 or near 1) -- a crude but conservative proxy for the
     probe not making contact or a bad reference-panel read.
  2. `qc_spectral_out_of_range()` -- flags physically-impossible
     reflectance values.
  3. `qc_spectral_group_outliers()` -- flags spectra whose *shape* (via
     spectral angle from the group median spectrum) is a robust outlier
     within their own group. Deliberately not PCA-based -- built for
     small replicate counts per group (e.g. 6), where per-group PCA is
     underpowered.

  **Known gap** (Henry, 2026-09-20): the most useful spectral check would
  be a moisture-index-style flag for incomplete drying, but validating one
  needs a dataset this hasn't had yet. Not attempted here.

- **Reporting/provenance** (`R/qc_report.R`): `qc_write_log()`,
  `qc_build_report()` (a per-check counts summary -- deliberately does
  *not* try to merge flags from different grains into one per-id column;
  see the function's own docstring for why that join has to live in the
  calling project), and `qc_apply_removals()` (the *only* function that
  mutates data -- everything else only flags; whether a flag becomes a
  removal is a per-project policy call, applied explicitly and
  separately).

## Per-project usage pattern

Nothing in this repo is deposit-ready by itself. A calling project:
1. Sources these files (or, once this is a real package, installs it).
2. Supplies its own config: which grouping variable applies to which
   trait, bounds tables, thresholds, and a remove-vs-flag policy per
   check -- committed with the project's own repo, since it's part of
   that deposit's provenance, not this library's.
3. Calls the `qc_*` functions per trait/check, writes results with
   `qc_write_log()`, and finishes with one `qc_build_report()` call.
4. Only if the project's policy calls for it: passes hard-bounds
   violations to `qc_apply_removals()`.

See the Mill Test 2024 repo (`Workflow_Calculations/code/`) for the first
worked example.

## Output contract

Every check returns (never mutates) a long-format tibble: an id column (or
columns) the caller supplied, the offending value(s), and a `check` /
`trait_flagged` column identifying what tripped. Nothing is dropped
silently -- flags accumulate into a report; only `qc_apply_removals()`
changes data, and only for checks a project has explicitly marked
"remove."
