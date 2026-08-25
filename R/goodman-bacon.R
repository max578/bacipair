#' Goodman-Bacon (2021) staggered-DiD decomposition
#'
#' Decomposes a two-way fixed-effects DiD estimator into a weighted
#' average of all 2x2 DiD estimators implied by the treatment timing,
#' as in Goodman-Bacon (2021). Wraps `bacondecomp::bacon()`; that
#' package is declared in `Suggests`, so the wrapper raises an
#' informative error if it is not installed.
#'
#' This decomposition is the diagnostic the translation memo (RDT-8)
#' relies on to flag "forbidden comparisons" - 2x2 cells in which an
#' already-treated cohort acts as the comparator for a not-yet-treated
#' cohort and so contaminates the TWFE coefficient with negative
#' weights. A clean decomposition (large weight on
#' `Treated vs Untreated`; small / non-negative weights on
#' `Earlier vs Later Treated`) is a precondition for trusting the
#' single-coefficient TWFE.
#'
#' Operates on the *long-format* treated-arm panel: only rows whose
#' `arm` column matches `treated_label` are passed to
#' `bacondecomp::bacon()`, so each observation is a (pair_id, year)
#' pair with a single 0/1 treatment indicator.
#'
#' @param panel A `data.table` containing the BACI panel (typically
#'   from [simulate_staggered_panel()] or a real-data port).
#' @param outcome Character. Outcome column name (default
#'   `"woody_cover"`).
#' @param id Character. Unit identifier column (default `"pair_id"`).
#' @param time Character. Time column (default `"year"`).
#' @param treated_label Character. Level of `arm` denoting treated
#'   units (default `"credited"`).
#' @param treatment_indicator Character. Name of the binary treatment
#'   indicator column (default `"post_treatment"`); must be 0/1 valued.
#' @return A list with components:
#'   \describe{
#'     \item{`decomposition`}{A `data.frame` returned by
#'       `bacondecomp::bacon()`; one row per 2x2 cell, with `treated`,
#'       `untreated`, `estimate`, `weight`, and `type`.}
#'     \item{`twfe_estimate`}{Numeric. The implied weighted-average
#'       TWFE estimate.}
#'     \item{`weight_summary`}{A `data.frame` summarising weights and
#'       weighted estimates by `type`.}
#'   }
#' @examples
#' if (requireNamespace("bacondecomp", quietly = TRUE)) {
#'   sim <- simulate_staggered_panel(seed = 20260520L, n_pairs = 30L,
#'                                   n_proponents = 8L)
#'   gb  <- goodman_bacon_decomp(sim$panel)
#'   gb$weight_summary
#' }
#' @export
goodman_bacon_decomp <- function(panel,
                                 outcome             = "woody_cover",
                                 id                  = "pair_id",
                                 time                = "year",
                                 treated_label       = "credited",
                                 treatment_indicator = "post_treatment") {
  stopifnot(
    data.table::is.data.table(panel),
    is.character(outcome), length(outcome) == 1L,
    is.character(id), is.character(time),
    is.character(treated_label),
    is.character(treatment_indicator),
    all(c("arm", id, time, outcome, treatment_indicator) %in% names(panel))
  )
  if (!requireNamespace("bacondecomp", quietly = TRUE)) {
    stop("goodman_bacon_decomp: package 'bacondecomp' is required ",
         "but not installed. Install it with ",
         "install.packages('bacondecomp').")
  }

  ## Local rebind to avoid colliding with the simulator's `treated_arm`
  ## column under data.table NSE (the column is a 0/1 integer, the
  ## function arg is a character level of the `arm` column).
  .keep_label <- treated_label
  d <- as.data.frame(panel[get("arm") == .keep_label,
                           c(id, time, outcome, treatment_indicator),
                           with = FALSE])
  if (!all(d[[treatment_indicator]] %in% c(0L, 1L, 0, 1))) {
    stop("goodman_bacon_decomp: '", treatment_indicator,
         "' must be 0/1 valued.")
  }

  fml <- stats::as.formula(
    paste(outcome, "~", treatment_indicator)
  )
  decomp <- bacondecomp::bacon(
    formula = fml,
    data    = d,
    id_var  = id,
    time_var = time,
    quietly = TRUE
  )

  ## bacondecomp::bacon() returns either a data.frame (no controls)
  ## with cols (treated, untreated, estimate, weight, type) or a list
  ## with $two_by_twos when controls are passed. We always call it
  ## without controls, so we expect the data.frame form.
  if (!is.data.frame(decomp)) {
    stop("goodman_bacon_decomp: unexpected return type from ",
         "bacondecomp::bacon().")
  }

  twfe_estimate <- sum(decomp$estimate * decomp$weight)
  weight_summary <- stats::aggregate(
    cbind(weight = decomp$weight,
          weighted_estimate = decomp$estimate * decomp$weight) ~
      type,
    data = decomp,
    FUN = sum
  )

  list(
    decomposition  = decomp,
    twfe_estimate  = twfe_estimate,
    weight_summary = weight_summary
  )
}
