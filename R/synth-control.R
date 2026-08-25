#' Project a vector onto the probability simplex
#'
#' Implements the Wang & Carreira-Perpinan (2013) projection: given a
#' real vector `v` of length `n`, returns the unique point in the
#' `n`-dimensional probability simplex (non-negative, sums to one) that
#' minimises Euclidean distance to `v`.
#'
#' @param v Numeric vector of length >= 1.
#' @return Numeric vector of the same length, non-negative, summing to 1.
#' @keywords internal
project_simplex <- function(v) {
  stopifnot(is.numeric(v), length(v) >= 1L, all(is.finite(v)))
  n <- length(v)
  u <- sort(v, decreasing = TRUE)
  cssv <- cumsum(u)
  rho_candidates <- which(u - (cssv - 1) / seq_len(n) > 0)
  if (length(rho_candidates) == 0L) {
    stop("project_simplex: no feasible rho found.")
  }
  rho <- rho_candidates[length(rho_candidates)]
  tau <- (cssv[rho] - 1) / rho
  pmax(v - tau, 0)
}

#' Solve a simplex-constrained least-squares problem
#'
#' Minimises `||X w - y||^2` over the probability simplex via projected
#' gradient with backtracking line search and opportunistic step growth.
#'
#' @param X Numeric matrix `[n x k]` of donor features.
#' @param y Numeric vector of length `n` of target values.
#' @param max_iter Integer (>= 10). Maximum iterations.
#' @param tol Numeric > 0. Convergence tolerance on objective improvement.
#' @param step_init Numeric > 0. Initial step size.
#' @return List with components `weights`, `objective`, `iterations`.
#' @keywords internal
solve_simplex_lsq <- function(X, y,
                              max_iter = 500L,
                              tol = 1e-9,
                              step_init = 1.0) {
  stopifnot(is.matrix(X), is.numeric(y),
            nrow(X) == length(y),
            ncol(X) >= 1L,
            max_iter >= 10L, tol > 0, step_init > 0,
            all(is.finite(X)), all(is.finite(y)))
  n_donors <- ncol(X)
  w <- rep(1 / n_donors, n_donors)

  XtX <- crossprod(X)
  Xty <- crossprod(X, y)

  obj <- function(w) sum((X %*% w - y) ^ 2)
  grad <- function(w) as.numeric(2 * (XtX %*% w - Xty))

  step <- step_init
  obj_prev <- obj(w)
  it <- 0L
  for (it in seq_len(max_iter)) {
    g <- grad(w)
    w_new <- project_simplex(w - step * g)
    obj_new <- obj(w_new)

    bt <- 0L
    while (obj_new > obj_prev + 1e-12 && bt < 30L) {
      step <- step * 0.5
      w_new <- project_simplex(w - step * g)
      obj_new <- obj(w_new)
      bt <- bt + 1L
    }

    delta <- obj_prev - obj_new
    w <- w_new
    obj_prev <- obj_new

    if (delta < tol) break
    if (bt == 0L) step <- step * 1.05
  }
  list(weights = w, objective = obj_prev, iterations = it)
}

#' Per-treated-unit synthetic-control weight fit
#'
#' Fits simplex-constrained donor weights for one treated unit, jointly
#' matching pre-treatment outcome trajectory and (standardised) static
#' predictors.
#'
#' @param treated_pre_outcomes Numeric vector indexed by pre-treatment year.
#' @param donor_pre_outcomes Numeric matrix `[years x donors]`.
#' @param treated_predictors Numeric vector of length `k` of static predictors.
#' @param donor_predictors Numeric matrix `[k x donors]`.
#' @param predictor_weight Numeric >= 0. Relative weight on the predictor
#'   block (default 1). Higher values force tighter predictor matching.
#' @return List with components `weights`, `objective`, `iterations`,
#'   `pre_residual`.
#' @keywords internal
fit_unit_synth_control <- function(treated_pre_outcomes,
                                   donor_pre_outcomes,
                                   treated_predictors,
                                   donor_predictors,
                                   predictor_weight = 1.0) {
  stopifnot(
    is.numeric(treated_pre_outcomes),
    is.matrix(donor_pre_outcomes),
    nrow(donor_pre_outcomes) == length(treated_pre_outcomes),
    is.numeric(treated_predictors),
    is.matrix(donor_predictors),
    nrow(donor_predictors) == length(treated_predictors),
    ncol(donor_pre_outcomes) == ncol(donor_predictors),
    predictor_weight >= 0
  )

  pred_sd <- apply(donor_predictors, 1L, stats::sd)
  pred_sd[pred_sd <= 0] <- 1
  donor_predictors_z <- (donor_predictors - rowMeans(donor_predictors)) /
    pred_sd
  treated_predictors_z <- (treated_predictors - rowMeans(donor_predictors)) /
    pred_sd

  X <- rbind(
    donor_pre_outcomes,
    sqrt(predictor_weight) * donor_predictors_z
  )
  y <- c(
    treated_pre_outcomes,
    sqrt(predictor_weight) * treated_predictors_z
  )

  fit <- solve_simplex_lsq(X, y)
  list(
    weights      = fit$weights,
    objective    = fit$objective,
    iterations   = fit$iterations,
    pre_residual = sqrt(mean(
      (donor_pre_outcomes %*% fit$weights -
         treated_pre_outcomes) ^ 2
    ))
  )
}

#' Fit a time-varying synthetic-control comparator on a BACI panel
#'
#' For each treated unit, constructs a per-unit synthetic-control
#' comparator by minimising the squared distance between the treated
#' pre-treatment trajectory + static predictors and a simplex-constrained
#' weighted average of every other pair's comparison-arm trajectory +
#' predictors. Returns per-(pair, year) ATT estimates and the
#' panel-average ATT.
#'
#' Pure base R; deterministic; no extra dependencies. The treated unit's
#' own matched 1:1 comparator is excluded from its donor pool to avoid
#' leakage; the pool is every *other* pair's comparison arm.
#'
#' @param panel A `data.table` with columns `pair_id`, `arm`, `year`,
#'   `treatment_year`, `post_treatment`, `woody_cover`, plus all entries
#'   of `predictor_cols`. An `ite` column (the simulator's oracle
#'   individual treatment effect) is read and reported back when
#'   present, but is never required -- the estimator runs on a real
#'   register that has no such column.
#' @param predictor_cols Character vector of predictor column names
#'   (default `c("baseline_cover", "elevation_m",
#'   "rainfall_climatology")`).
#' @param predictor_weight Numeric >= 0. Relative weight on the
#'   predictor block (default 20).
#' @param min_pre_years Integer >= 2. Minimum pre-treatment years
#'   required per treated unit (default 3).
#' @return List with components:
#'   \describe{
#'     \item{`fits`}{List of per-unit fits.}
#'     \item{`per_unit_att`}{`data.table` of per-unit estimated ATT,
#'       plus the oracle ATT (`oracle_unit`) when the input panel
#'       carried an `ite` column.}
#'     \item{`panel_att_hat`}{Numeric. Panel-average estimated ATT.}
#'     \item{`pre_fit_rmse`}{Numeric. Root-mean-square pre-fit residual
#'       across treated units.}
#'   }
#' @examples
#' sim <- simulate_staggered_panel(seed = 20260520L, n_pairs = 30L,
#'                                  n_proponents = 8L)
#' fit <- synth_control_baci(sim$panel)
#' fit$panel_att_hat
#' @export
synth_control_baci <- function(panel,
                               predictor_cols = c("baseline_cover",
                                                  "elevation_m",
                                                  "rainfall_climatology"),
                               predictor_weight = 20.0,
                               min_pre_years = 3L) {
  stopifnot(
    data.table::is.data.table(panel),
    all(c("pair_id", "year", "arm", "treatment_year",
          "post_treatment", "woody_cover",
          predictor_cols) %in% names(panel)),
    min_pre_years >= 2L,
    predictor_weight >= 0
  )

  treated_pairs <- unique(panel[get("arm") == "credited", get("pair_id")])
  comparison_long <- panel[get("arm") == "comparison"]
  donors <- unique(comparison_long$pair_id)

  if (length(donors) < 5L) {
    abstain(
      paste0("synth_control_baci: donor pool has <5 comparison pairs ",
             "(", length(donors), "); the simplex-constrained fit is not ",
             "identified below this floor."),
      reason = "donor_pool_too_small",
      data   = list(n_donors = length(donors), floor = 5L)
    )
  }

  predictors_dt <- unique(panel[, c("pair_id", predictor_cols),
                                with = FALSE])

  donor_predictors <- t(as.matrix(
    predictors_dt[match(donors, get("pair_id")),
                  predictor_cols, with = FALSE]
  ))

  donor_panel_wide <- data.table::dcast(
    comparison_long, year ~ pair_id,
    value.var = "woody_cover", fun.aggregate = mean
  )
  donor_years <- donor_panel_wide$year
  donor_outcome_mat <- as.matrix(
    donor_panel_wide[, -"year", with = FALSE]
  )
  donor_outcome_mat <- donor_outcome_mat[, as.character(donors),
                                         drop = FALSE]

  results <- vector("list", length(treated_pairs))
  for (k in seq_along(treated_pairs)) {
    i <- treated_pairs[[k]]
    t_g <- panel[get("pair_id") == i & get("arm") == "credited",
                 unique(get("treatment_year"))]
    if (length(t_g) != 1L) {
      stop("synth_control_baci: pair ", i,
           " has non-unique treatment_year.")
    }
    pre_mask <- donor_years < t_g
    if (sum(pre_mask) < min_pre_years) {
      abstain(
        paste0("synth_control_baci: pair ", i, " has only ",
               sum(pre_mask), " pre-treatment years (< min_pre_years = ",
               min_pre_years, "); the pre-treatment fit is not ",
               "identified below this floor."),
        reason = "insufficient_pre_years",
        data   = list(pair_id = i, n_pre_years = sum(pre_mask),
                      min_pre_years = min_pre_years)
      )
    }

    treated_outcomes <- panel[
      get("pair_id") == i & get("arm") == "credited"
    ][order(get("year"))]
    treated_pre <- treated_outcomes[get("year") < t_g, get("woody_cover")]

    treated_predictors <- as.numeric(
      predictors_dt[get("pair_id") == i, predictor_cols, with = FALSE]
    )

    donor_pre <- donor_outcome_mat[pre_mask, , drop = FALSE]

    keep_donors <- donors != i
    donor_pre <- donor_pre[, keep_donors, drop = FALSE]
    donor_predictors_i <- donor_predictors[, keep_donors, drop = FALSE]
    donor_full_i <- donor_outcome_mat[, keep_donors, drop = FALSE]

    fit <- fit_unit_synth_control(
      treated_pre_outcomes = treated_pre,
      donor_pre_outcomes   = donor_pre,
      treated_predictors   = treated_predictors,
      donor_predictors     = donor_predictors_i,
      predictor_weight     = predictor_weight
    )

    counterfactual <- as.numeric(donor_full_i %*% fit$weights)
    counterfactual_dt <- data.table::data.table(
      pair_id        = i,
      year           = donor_years,
      counterfactual = counterfactual
    )
    ## `ite` is the simulator's oracle-only column (present only on a
    ## simulated panel, never on a real register) -- attach it when
    ## present so oracle grading still works on simulator output, but
    ## never require it: the estimated ATT is defined purely from
    ## `woody_cover - counterfactual` and must run on real data (F1).
    has_oracle <- "ite" %in% names(treated_outcomes)
    keep_cols <- c("pair_id", "year", "woody_cover",
                   "treatment_year", "post_treatment")
    if (has_oracle) {
      keep_cols <- c(keep_cols, "ite")
    }
    fit_dt <- merge(
      treated_outcomes[, keep_cols, with = FALSE],
      counterfactual_dt,
      by = c("pair_id", "year"),
      all.x = TRUE
    )
    fit_dt[, "att_hat" := get("woody_cover") - get("counterfactual")]

    results[[k]] <- list(
      pair_id        = i,
      treatment_year = t_g,
      weights        = fit$weights,
      pre_residual   = fit$pre_residual,
      iterations     = fit$iterations,
      panel          = fit_dt
    )
  }

  per_unit <- data.table::rbindlist(lapply(results, `[[`, "panel"),
                                    fill = TRUE)
  panel_att <- per_unit[get("post_treatment") == 1L, mean(get("att_hat"))]
  if ("ite" %in% names(per_unit)) {
    per_unit_att <- per_unit[get("post_treatment") == 1L,
                             list(att_hat_unit = mean(get("att_hat")),
                                  oracle_unit  = mean(get("ite"))),
                             by = "pair_id"]
  } else {
    per_unit_att <- per_unit[get("post_treatment") == 1L,
                             list(att_hat_unit = mean(get("att_hat"))),
                             by = "pair_id"]
  }

  list(
    fits           = results,
    per_unit_att   = per_unit_att,
    panel_att_hat  = panel_att,
    pre_fit_rmse   = sqrt(mean(
      vapply(results, function(r) r$pre_residual ^ 2, numeric(1))
    ))
  )
}
