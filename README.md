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

> **Physical bounds and literature bounds are not the same thing**
> (Henry, 2026-09-23). Only *physically impossible* values — %C outside
> 0–100, a water-content fraction ≥ 1 — are candidates for removal. A
> literature or TRY range is **flag-only, always**: this lab works with very
> diverse material, and a true value outside TRY's compiled range is not
> surprising. Keep the two in separate tables in the project config, and
> never put a literature-sourced trait in `removal_policy`.

- **Trait QC** (`R/trait_qc.R`):
  1. `qc_hard_bounds()` -- flag values outside a bound. Used for both bound
     types; what differs is the project's policy on what happens next (see
     the note above). Catches mis-entered data: decimal shifts, sign errors,
     unit mismatches.
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
  6. `qc_relative_error()` -- flags a result whose own reported dispersion
     (an SD over the assay's duplicate/triplicate measures) is too large a
     fraction of the value. Asks "did the assay agree with itself", which no
     comparison against other samples can answer.
  7. `qc_blank_drift()` -- flags a QC blank that **gains** weight across a
     sequence of gravimetric processing steps. A blank should only ever lose
     weight; a gain means it leaked and took on material, so it was shedding
     sample material too. This is a **batch-level** signal -- the caller is
     responsible for propagating it to everything that shared the batch.
  8. `qc_sequential_cascade()` -- propagates a failure forward through an
     assay whose stages each operate on the residue of the last (ANKOM fiber:
     NDF, then ADF, then ADL). Deliberately precautionary: it flags on
     position in the chain, not on evidence the downstream value is itself
     wrong, so check the assay's arithmetic before removing anything on this
     basis.
  9. `qc_derived_mismatch()` -- compares a spreadsheet-reported derived value
     against the same quantity recomputed in code from the raw measurements
     it comes from. **This catches a failure mode nothing else here can see.**
     When results arrive as a workbook whose derived columns are formulas over
     typed-in measurements, the arithmetic itself can break while the output
     stays perfectly plausible: a dragged fill-handle or an inserted row
     shifts a relative cell reference, so a formula silently reads its
     neighbour's weight. The result sits inside every bound and matches its
     own batch, so bounds, group, distribution and literature checks all pass
     it. The Mill Test 2024 ANKOM sheets carried exactly this -- one corrupted
     value was impossible and got caught, a second was an 11.6 mg/g lignin
     error that looked entirely normal. It also catches a stale cached value,
     where the formula is right but the file was saved without recalculating
     (readers like `readxl` return the cache, not the formula). Worth running
     wherever raw inputs sit alongside derived outputs. The stronger move,
     where the protocol is documented, is to **carry the recomputed value
     forward** and let this check exist only to report the sheet's defects --
     then a future reference slip cannot change a published number.
  10. `qc_try_bounds()` (`R/try_reference.R`) -- pulls min/max literature
     ranges for a set of traits out of a bulk TRY database export, to be
     used as **flag-only** bounds. Never feed these into a removal policy
     (see the note above). Check `UnitName` on what comes back: some TRY
     `TraitName`s pool incompatible units in `StdValue`, which produces a
     meaningless range if you take min/max blindly.

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
     spectral angle) is a robust outlier within their own group. The
     reference is always a **real measured scan**, picked either as the
     `medoid` (smallest summed angle to the others — best for controlled
     lab settings) or by `albedo` (the scan closest to the group's median
     brightness — expected to matter more for fresh-leaf and field work).
     Set `reference` per project; the method and the scan used are both
     returned so they land in the QC log. Deliberately not PCA-based —
     built for small replicate counts per group, where per-group PCA is
     underpowered.
  4. `qc_spectral_band_threshold()` -- flags spectra whose reflectance at
     one band (or averaged over a band range) is implausible. The SHIFT
     dried-ground DAAC product drops everything with R350 ≥ 0.9; dry leaf
     powder sits far below that, so a value that high means the white
     reference or an empty puck got measured. More direct than
     `qc_spectral_flat()` for that specific failure — keep both.
  5. `qc_spectral_pca_flags()` -- flags spectra **for manual review** at two
     levels: each replicate against its own sample centroid, and each sample
     centroid against its group's. The PCA is fit **once over the whole
     dataset** and both levels are measured in that shared space, which is
     what makes it work on thin replicate sets. Distances use the top `n_pc`
     components scaled by their global MAD; z-scores are pooled across all
     replicates/samples, because a single group is far too small to
     robust-z within. Review-only — never wire this to a removal policy.
  6. `qc_spectral_splice_jump()` -- flags unusually large steps at detector
     boundaries, measured **before** jump correction (afterwards the
     evidence is gone by construction). Matters most for field spectra,
     where detector temperature moves the jumps. `qc_detect_jump_corrected()`
     tells you whether a file has already been corrected — a
     boundary-matching correction leaves `R(splice) == R(splice+1)` exactly,
     which is otherwise vanishingly rare. **A file named "raw" is not
     necessarily uncorrected; check before trusting it.**
     `qc_splice_defaults()` holds ASD and SVC/.sig boundaries and treats
     `"svc"` and `"sig"` as the same instrument.

  **Known gap** (Henry, 2026-09-20): the most useful spectral check would
  be a moisture-index-style flag for incomplete drying, but validating one
  needs a dataset this hasn't had yet. An experiment to produce that
  dataset is designed (paired oven trays, one scanned and one weighed) but
  not yet run — see the Mill Test CLAUDE.md. Not attempted here.

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
