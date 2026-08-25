#' bacipair: BACI / paired-DiD / synthetic-control for carbon-MRV audits
#'
#' A reproducibility-by-construction toolkit for Before-After-Control-Impact
#' (BACI) audits of carbon-MRV registers. Designed against the Macintosh et al.
#' (2024) Australian Human-Induced Regeneration setting but applicable to any
#' paired BACI register where treated and matched comparator units are observed
#' over a panel of years.
#'
#' @section Public API:
#'
#' \describe{
#'   \item{[simulate_staggered_panel()]}{Staggered-treatment HIR-style BACI
#'     panel simulator with explicit oracle panel-average and per-unit ATT.}
#'   \item{[baci_paired_did()]}{Paired difference-in-differences on a BACI
#'     panel with proponent-cluster-robust standard errors (CR1).}
#'   \item{[synth_control_baci()]}{Time-varying synthetic-control comparator
#'     with simplex-constrained donor weights (Wang and Carreira-Perpinan,
#'     2013) and rainfall / elevation / baseline-cover predictors.}
#'   \item{[proponent_cluster_bootstrap()]}{Proponent-level resampling
#'     bootstrap for the synthetic-control or paired-DiD point estimate.}
#'   \item{[goodman_bacon_decomp()]}{Goodman-Bacon (2021) two-by-two
#'     decomposition of staggered-DiD weights, exposed as a guarded wrapper
#'     around `bacondecomp::bacon()`.}
#'   \item{[bacipair_targets_plan()]}{`targets`-skeleton helper that returns
#'     a list of `targets::tar_target` calls covering simulate -> SC fit ->
#'     paired-DiD fit -> bootstrap -> Goodman-Bacon decomposition ->
#'     oracle-recovery JSON log.}
#' }
#'
#' @section Reproducibility:
#'
#' Every public function is deterministic given a seed. No silent error
#' handling: invalid inputs trigger `stop()` with diagnostic messages.
#'
#' @section References:
#'
#' Goodman-Bacon, A. (2021). Difference-in-differences with variation in
#' treatment timing. *Journal of Econometrics*, 225(2), 254-277.
#' \doi{10.1016/j.jeconom.2021.03.014}.
#'
#' Macintosh, A., Butler, D., Larraondo, P., Evans, M. C., Ansell, D.,
#' Waschka, M., Fensham, R., Eldridge, D., Lindenmayer, D., Gibbons, P., &
#' Summerfield, P. (2024). Australian human-induced native forest
#' regeneration carbon offset projects have limited impact on changes in
#' woody vegetation cover and carbon removals. *Communications Earth &
#' Environment*, 5, Article 149. \doi{10.1038/s43247-024-01313-x}.
#'
#' Wang, W., & Carreira-Perpinan, M. A. (2013). Projection onto the
#' probability simplex: An efficient algorithm with a simple proof, and an
#' application. *arXiv preprint*. \doi{10.48550/arXiv.1309.1541}.
#'
#' @keywords internal
#' @importFrom data.table data.table is.data.table copy CJ rbindlist dcast
#'   setcolorder :=  .N
#' @importFrom stats aggregate as.formula qnorm rbinom rnorm sd quantile
#' @importFrom utils packageVersion sessionInfo capture.output globalVariables
"_PACKAGE"

## Suppress R CMD check NOTEs about data.table NSE column names used inside
## j-expressions. These are referenced by character-string column names from
## function arguments and are fully validated upstream.
utils::globalVariables(c(
  "pair_id", "year", "arm", "treatment_year", "post_treatment", "treated_arm",
  "post_idx", "ite", "noise", "woody_cover", "site_pair_shift",
  "rainfall_anomaly", "att_hat", "counterfactual",
  "true_additionality_max", "has_real_additionality",
  "baseline_cover", "elevation_m", "rainfall_climatology",
  "rainfall_sensitivity", "proponent_id"
))
