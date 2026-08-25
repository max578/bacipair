#' Proponent-level cluster bootstrap
#'
#' Resamples whole proponent clusters (with replacement) from a BACI panel
#' and refits a chosen estimator on each draw. Returns the bootstrap
#' distribution of the panel-average ATT, along with a 95% percentile
#' interval and the bootstrap standard error.
#'
#' Cluster-level resampling is the right inferential unit when one
#' proponent runs many adjacent projects (Macintosh et al., 2024) -
#' resampling pairs would understate dependence.
#'
#' @param panel A `data.table` with the columns required by the chosen
#'   estimator, plus the cluster column (default `proponent_id`).
#' @param estimator Character, one of `"synth_control"` (default) or
#'   `"paired_did"`. Selects which estimator to refit on each bootstrap
#'   draw.
#' @param n_boot Integer (>= 50). Number of bootstrap draws (default 200).
#' @param seed Integer. RNG seed (default `20260521L`).
#' @param cluster_col Character. Cluster column name (default
#'   `"proponent_id"`).
#' @param ... Extra arguments forwarded to the chosen estimator. For
#'   `"synth_control"`, see [synth_control_baci()] (`predictor_cols`,
#'   `predictor_weight`, `min_pre_years`). For `"paired_did"`, see
#'   [baci_paired_did()] (`outcome`, `treated_label`, `comparison_label`).
#' @return List with components:
#'   \describe{
#'     \item{`draws`}{Numeric vector of length `n_boot` of point
#'       estimates per draw.}
#'     \item{`se`}{Numeric. Bootstrap standard deviation across draws.}
#'     \item{`ci_lo`,`ci_hi`}{Numeric. 2.5th / 97.5th percentile.}
#'     \item{`n_boot`,`estimator`,`seed`}{Bookkeeping.}
#'   }
#' @examples
#' sim <- simulate_staggered_panel(seed = 20260520L, n_pairs = 30L,
#'                                  n_proponents = 8L)
#' boot <- proponent_cluster_bootstrap(sim$panel, estimator = "paired_did",
#'                                     n_boot = 50L)
#' boot$ci_lo
#' boot$ci_hi
#' @export
proponent_cluster_bootstrap <- function(panel,
                                        estimator   = c("synth_control",
                                                        "paired_did"),
                                        n_boot      = 200L,
                                        seed        = 20260521L,
                                        cluster_col = "proponent_id",
                                        ...) {
  estimator <- match.arg(estimator)
  stopifnot(
    data.table::is.data.table(panel),
    cluster_col %in% names(panel),
    n_boot >= 50L,
    is.numeric(seed), length(seed) == 1L
  )

  set.seed(seed)
  clusters <- unique(panel[[cluster_col]])
  draws <- numeric(n_boot)

  for (b in seq_len(n_boot)) {
    sampled <- sample(clusters, size = length(clusters), replace = TRUE)
    boot_panel <- data.table::rbindlist(
      lapply(seq_along(sampled), function(idx) {
        sub <- panel[get(cluster_col) == sampled[idx]]
        sub <- data.table::copy(sub)
        sub[, "pair_id" := get("pair_id") + idx * 1e6L]
        sub
      })
    )
    draws[b] <- if (estimator == "synth_control") {
      synth_control_baci(boot_panel, ...)$panel_att_hat
    } else {
      baci_paired_did(boot_panel, ...)$att_hat
    }
  }

  list(
    draws     = draws,
    se        = stats::sd(draws),
    ci_lo     = stats::quantile(draws, 0.025, names = FALSE),
    ci_hi     = stats::quantile(draws, 0.975, names = FALSE),
    n_boot    = n_boot,
    estimator = estimator,
    seed      = seed
  )
}
