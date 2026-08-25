# bacipair (development version)

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

* Package is now under git version control (`max578/bacipair`, private)
  and carries a standard `R-CMD-check` GitHub Actions workflow (F10).

# bacipair 0.1.0 (2026-05-03)

Initial scaffold release.

## Public API

* `simulate_staggered_panel()` — staggered-treatment HIR-style BACI
  panel simulator with explicit oracle panel-average and per-unit ATT.
* `baci_paired_did()` — paired difference-in-differences on a BACI
  panel with proponent-cluster-robust standard errors (CR1).
* `synth_control_baci()` — time-varying synthetic-control comparator
  with simplex-constrained donor weights and rainfall / elevation /
  baseline-cover predictors.
* `proponent_cluster_bootstrap()` — proponent-level resampling
  bootstrap for the synthetic-control or paired-DiD point estimate.
* `goodman_bacon_decomp()` — Goodman-Bacon (2021) two-by-two
  decomposition of staggered-DiD weights, exposed as a guarded
  wrapper around `bacondecomp::bacon()`.
* `bacipair_targets_plan()` — `targets`-skeleton helper that returns
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
