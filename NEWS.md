# bacipair (development version)

## Documentation and CI

* The package has a logo.
* README carries status badges and lists related packages.
* CI runs lintr (report-only) and test coverage.

## Bug fixes

* `synth_control_baci()` no longer requires the simulator-only oracle
  column `ite`; it runs on a real panel and attaches `ite`/`oracle_unit`
  only when the input carries them (F1).

## Contract

* The synthetic-control estimator's donor-pool-too-small and
  insufficient-pre-treatment-years declines are now a classed
  `bacipair_abstention` condition (also carrying `orchestra_refusal`),
  per the orchestra's refusal/abstention naming convention, rather than
  a bare `stop()` (F4).

## Infrastructure

* Package is now under git version control and carries a standard `R-CMD-check` GitHub Actions workflow (F10).

## Documentation

* Added a README with a three-line purpose statement, install
  instructions and a minimal runnable example (F8).
* Added an "Estimand and estimators" section to the getting-started
  vignette giving the ATT definition, the paired-DiD two-way
  fixed-effects specification, the CR1 cluster-robust variance formula
  and the simplex-constrained synthetic-control objective, all
  previously stated in prose only with no checkable mathematics (F9).
* Added three captioned, interpreted figures to the getting-started
  vignette: pre-treatment fit and post-treatment gap for the best- and
  worst-fitting treated units, per-unit estimated ATT against the
  simulator oracle, and the donor-weight matrix (F12).
* Restated the vignette's kill-gate claim: it now says plainly that
  the shipped `synth_control_baci()` is the uncorrected estimator, and
  quotes the uncorrected audit-grade bias (0.0347 pp against the
  0.15 pp floor) rather than implying reproduction of a corrected
  extension that is not in the package API (F2).
* Inlined the kill-gate floor and achieved-bias numbers directly in
  the vignette text instead of pointing the reader at
  `deliverables/RDT-20/...` paths that do not ship with the package
  (F13).
* Replaced the vignette's unsupported "small but real residual bias"
  characterisation of the paired-DiD estimator with the measured
  relative bias (37 per cent of the oracle effect, about twice the
  cluster-robust standard error) (F14).
* Expanded the vignette's unexplained `CER`/`ACCU`/`HIR`/`MRV`
  acronym cluster on first use, and removed a `+` used as a prose
  conjunction.
* Corrected `@param` documentation of `additionality_mean_pp` /
  `additionality_sd_pp` in `simulate_staggered_panel()`: the draw is
  an untruncated Normal, not truncated-Normal as previously
  documented, and can realise a negative additionality for a pair
  flagged as genuinely additional (F3, documentation only -- the
  simulator's numerical behaviour is unchanged in this release).
* Corrected the `@return` documentation of `baci_paired_did()`'s
  `coefficients` element: the returned vector holds the ATT
  interaction and non-base-year dummies on the demeaned design, not
  an intercept, treated-arm dummy or post-treatment dummy, all three
  of which are absorbed by demeaning and were never returned (F16).
* Added `@examples` to the two exported functions that previously had
  none, `is_bacipair_abstention()` and `write_oracle_recovery_log()`.

# bacipair 0.1.0 (2026-05-03)

Initial scaffold release.

## Public API

* `simulate_staggered_panel()` -- staggered-treatment HIR-style BACI
  panel simulator with explicit oracle panel-average and per-unit ATT.
* `baci_paired_did()` -- paired difference-in-differences on a BACI
  panel with proponent-cluster-robust standard errors (CR1).
* `synth_control_baci()` -- time-varying synthetic-control comparator
  with simplex-constrained donor weights and rainfall / elevation /
  baseline-cover predictors.
* `proponent_cluster_bootstrap()` -- proponent-level resampling
  bootstrap for the synthetic-control or paired-DiD point estimate.
* `goodman_bacon_decomp()` -- Goodman-Bacon (2021) two-by-two
  decomposition of staggered-DiD weights, exposed as a guarded
  wrapper around `bacondecomp::bacon()`.
* `bacipair_targets_plan()` -- `targets`-skeleton helper that returns
  a list of `targets::tar_target` calls covering simulate -> SC fit ->
  paired-DiD fit -> bootstrap -> Goodman-Bacon decomposition ->
  oracle-recovery JSON log.

## Reproducibility

* Every public function is deterministic given a seed.
* No silent error handling: invalid inputs trigger `stop()` with
  diagnostic messages.
* `R CMD check --as-cran` clean (zero warnings, zero notes other
  than first-release notes).
* Vignette `getting-started.Rmd` covers the synthetic
  ground-truth-recovery audit and the published kill-gate JSON log.

## Citation

If you use this package, please cite both the package and the source
methodology paper:

> Macintosh, A., et al. (2024). Australian human-induced native
> forest regeneration carbon offset projects have limited impact on
> changes in woody vegetation cover and carbon removals.
> *Communications Earth & Environment*, 5, Article 149.
> <https://doi.org/10.1038/s43247-024-01313-x>
