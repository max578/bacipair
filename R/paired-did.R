#' Paired difference-in-differences on a BACI panel
#'
#' Estimates the average treatment effect on the treated (ATT) on a paired
#' BACI panel via a two-way fixed-effects regression of the outcome on
#' (a) a treated-arm-by-post-treatment interaction, plus (b) pair fixed
#' effects via within-pair-and-arm demeaning, and (c) year fixed effects.
#'
#' Standard errors are cluster-robust at the proponent level (default) or
#' at the pair level if `cluster_col = "pair_id"`. The estimator is the
#' canonical CR1 sandwich form used in econometric panel work
#' (Cameron & Miller, 2015), and is computed without external
#' dependencies.
#'
#' Operationally, this is the *naive* paired-DiD comparator the
#' translation memo (RDT-8) uses as the foil for the synthetic-control
#' estimator: when elevation / rainfall climatology covary with treatment
#' timing, paired-DiD has residual bias that the SC estimator nets out.
#'
#' @param panel A `data.table` containing one row per (pair, year, arm)
#'   observation. Must contain `pair_id`, `arm`, `year`, `treatment_year`,
#'   and the outcome / cluster columns.
#' @param outcome Character. Outcome column name (default `"woody_cover"`).
#' @param treated_label,comparison_label Character scalars. Levels of
#'   the `arm` column denoting treated vs comparison units (defaults
#'   `"credited"` / `"comparison"`).
#' @param cluster_col Character. Cluster identifier (default
#'   `"proponent_id"`); must be present in the panel.
#' @return List with components:
#'   \describe{
#'     \item{`att_hat`}{Numeric. Estimated ATT (paired-DiD coefficient).}
#'     \item{`se_cluster`}{Numeric. Cluster-robust (CR1) standard error.}
#'     \item{`t_stat`}{Numeric. Robust t-statistic.}
#'     \item{`ci_lo`,`ci_hi`}{Numeric. 95% normal-approximation CI.}
#'     \item{`n_obs`,`n_clusters`,`n_pairs`,`n_years`}{Integers. Sample
#'       counts.}
#'     \item{`coefficients`}{Numeric named vector of all estimated
#'       coefficients (intercept, treated-arm dummy, post-treatment
#'       dummy, ATT interaction, year FE).}
#'     \item{`call`}{Captured `match.call()`.}
#'   }
#' @examples
#' sim <- simulate_staggered_panel(seed = 20260520L, n_pairs = 30L,
#'                                  n_proponents = 8L)
#' did <- baci_paired_did(sim$panel)
#' did$att_hat
#' did$se_cluster
#' @export
baci_paired_did <- function(panel,
                            outcome          = "woody_cover",
                            treated_label    = "credited",
                            comparison_label = "comparison",
                            cluster_col      = "proponent_id") {
  stopifnot(
    data.table::is.data.table(panel),
    is.character(outcome), length(outcome) == 1L,
    is.character(treated_label), is.character(comparison_label),
    is.character(cluster_col), length(cluster_col) == 1L,
    all(c("pair_id", "arm", "year", "treatment_year",
          outcome, cluster_col) %in% names(panel))
  )

  arm_vals <- unique(panel$arm)
  if (!all(c(treated_label, comparison_label) %in% arm_vals)) {
    stop("baci_paired_did: panel$arm must contain both '", treated_label,
         "' and '", comparison_label, "'.")
  }

  ## Pre-extract arm levels and re-bind to local names that don't collide
  ## with simulator-output column names (e.g. `treated_arm`) under
  ## data.table NSE.
  .keep_arms <- c(treated_label, comparison_label)
  d <- data.table::copy(panel[get("arm") %in% .keep_arms])
  d[, "treated_arm_i" := as.integer(get("arm") == treated_label)]
  d[, "post_i" := as.integer(get("year") >= get("treatment_year"))]
  d[, "did_i" := get("treated_arm_i") * get("post_i")]

  ## Standard TWFE staggered-DiD specification: pair-and-arm fixed
  ## effects (absorbed by within-group demeaning), year fixed effects
  ## (expanded as dummies, base year dropped), and a single
  ## treated-arm-by-post-treatment interaction whose coefficient is the
  ## ATT. We omit `treated_arm` and `post` as separate regressors
  ## because they are absorbed by the pair-arm and year FE blocks
  ## respectively under this design.

  years_obs <- sort(unique(d$year))
  if (length(years_obs) < 2L) {
    stop("baci_paired_did: panel must contain >= 2 distinct years.")
  }
  year_levels <- years_obs[-1L]  # drop first as base

  for (y in year_levels) {
    col_y <- paste0("yr_", y)
    d[, (col_y) := as.integer(get("year") == y)]
  }

  reg_cols <- c(outcome, "did_i", paste0("yr_", year_levels))

  for (col in reg_cols) {
    d[, (paste0(col, "_dm")) :=
        get(col) - mean(get(col)),
      by = c("pair_id", "arm")]
  }

  X_dm <- as.matrix(d[, paste0(c("did_i", paste0("yr_", year_levels)),
                                "_dm"),
                       with = FALSE])
  y_dm <- as.matrix(d[, paste0(outcome, "_dm"), with = FALSE])

  XtX <- crossprod(X_dm)
  Xty <- crossprod(X_dm, y_dm)

  qr_obj <- qr(XtX)
  if (qr_obj$rank < ncol(XtX)) {
    stop("baci_paired_did: design matrix is rank-deficient (rank ",
         qr_obj$rank, " of ", ncol(XtX),
         "). Check year coverage / missing data.")
  }
  beta <- as.numeric(solve(qr_obj, Xty))
  names(beta) <- colnames(X_dm)

  resid <- as.numeric(y_dm - X_dm %*% beta)
  d[, "_resid" := resid]

  ## CR1 cluster-robust variance: V = bread %*% meat %*% bread,
  ## bread = (X'X)^-1, meat = sum_g X_g' u_g u_g' X_g, with the standard
  ## degrees-of-freedom adjustment c = G/(G-1) * (N-1)/(N-K).
  G <- length(unique(d[[cluster_col]]))
  N <- nrow(X_dm)
  K <- ncol(X_dm)
  c_adj <- (G / (G - 1)) * ((N - 1) / (N - K))

  bread <- solve(XtX)
  meat <- matrix(0, nrow = K, ncol = K)
  for (g in unique(d[[cluster_col]])) {
    idx <- which(d[[cluster_col]] == g)
    Xg <- X_dm[idx, , drop = FALSE]
    ug <- resid[idx]
    Xg_ug <- crossprod(Xg, ug)
    meat <- meat + Xg_ug %*% t(Xg_ug)
  }
  V <- c_adj * (bread %*% meat %*% bread)

  did_idx <- which(colnames(X_dm) == "did_i_dm")
  att_hat <- beta[[did_idx]]
  se_cluster <- sqrt(V[did_idx, did_idx])
  t_stat <- att_hat / se_cluster
  z975 <- stats::qnorm(0.975)

  list(
    att_hat      = att_hat,
    se_cluster   = se_cluster,
    t_stat       = t_stat,
    ci_lo        = att_hat - z975 * se_cluster,
    ci_hi        = att_hat + z975 * se_cluster,
    n_obs        = N,
    n_clusters   = G,
    n_pairs      = length(unique(d$pair_id)),
    n_years      = length(years_obs),
    coefficients = beta,
    call         = match.call()
  )
}
