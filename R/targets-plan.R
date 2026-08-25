#' Build a `targets`-skeleton plan for a bacipair audit
#'
#' Returns a list of `targets::tar_target` calls covering the canonical
#' pipeline:
#'
#' 1. simulate the staggered-treatment panel (or load a real BACI panel),
#' 2. fit the synthetic-control estimator,
#' 3. fit the paired-DiD estimator,
#' 4. proponent-cluster bootstrap of the SC point estimate,
#' 5. Goodman-Bacon decomposition diagnostic,
#' 6. write a JSON oracle-recovery log.
#'
#' Skeleton means: this is a list ready to drop into `targets::tar_plan()`
#' (or `targets::tar_pipeline()`). The function does *not* execute the
#' pipeline. To run it, paste the returned list into a `_targets.R`
#' script and call `targets::tar_make()`.
#'
#' The `targets` package is in `Suggests`, so this function raises an
#' informative error if `targets` is not installed.
#'
#' @param simulate_args A named list passed to
#'   [simulate_staggered_panel()] (default `list(seed = 20260520L,
#'   n_pairs = 200L)`).
#' @param synth_args A named list passed to [synth_control_baci()]
#'   (default empty).
#' @param did_args A named list passed to [baci_paired_did()] (default
#'   empty).
#' @param boot_args A named list passed to
#'   [proponent_cluster_bootstrap()] (default `list(n_boot = 200L,
#'   seed = 20260521L)`).
#' @param log_path Character. Path for the JSON oracle-recovery log
#'   (default `"logs/oracle-recovery.json"`).
#' @return A named list of `targets::tar_target()` calls, one per
#'   pipeline node. Names are `panel`, `synth_control_fit`,
#'   `paired_did_fit`, `cluster_boot`, `goodman_bacon`,
#'   `oracle_recovery_log`.
#' @examples
#' if (requireNamespace("targets", quietly = TRUE)) {
#'   plan <- bacipair_targets_plan()
#'   names(plan)
#' }
#' @export
bacipair_targets_plan <- function(simulate_args = list(seed = 20260520L,
                                                       n_pairs = 200L),
                                  synth_args    = list(),
                                  did_args      = list(),
                                  boot_args     = list(n_boot = 200L,
                                                       seed = 20260521L),
                                  log_path      = "logs/oracle-recovery.json") {
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("bacipair_targets_plan: package 'targets' is required ",
         "but not installed. Install it with ",
         "install.packages('targets').")
  }
  stopifnot(
    is.list(simulate_args), is.list(synth_args),
    is.list(did_args), is.list(boot_args),
    is.character(log_path), length(log_path) == 1L
  )

  list(
    panel = targets::tar_target_raw(
      "panel",
      substitute(do.call(bacipair::simulate_staggered_panel, args),
                 list(args = simulate_args))
    ),
    synth_control_fit = targets::tar_target_raw(
      "synth_control_fit",
      substitute(do.call(bacipair::synth_control_baci,
                         c(list(panel = panel$panel), args)),
                 list(args = synth_args))
    ),
    paired_did_fit = targets::tar_target_raw(
      "paired_did_fit",
      substitute(do.call(bacipair::baci_paired_did,
                         c(list(panel = panel$panel), args)),
                 list(args = did_args))
    ),
    cluster_boot = targets::tar_target_raw(
      "cluster_boot",
      substitute(do.call(bacipair::proponent_cluster_bootstrap,
                         c(list(panel = panel$panel,
                                estimator = "synth_control"),
                           args)),
                 list(args = boot_args))
    ),
    goodman_bacon = targets::tar_target_raw(
      "goodman_bacon",
      quote(bacipair::goodman_bacon_decomp(panel$panel))
    ),
    oracle_recovery_log = targets::tar_target_raw(
      "oracle_recovery_log",
      substitute(bacipair::write_oracle_recovery_log(
        panel = panel,
        synth_control_fit = synth_control_fit,
        paired_did_fit = paired_did_fit,
        cluster_boot = cluster_boot,
        path = path
      ), list(path = log_path))
    )
  )
}

#' Write the oracle-recovery JSON log
#'
#' Captures the comparison between the simulated oracle ATT and the
#' point estimates from the synthetic-control and paired-DiD estimators
#' (plus the cluster-bootstrap CI on the SC estimator), in the schema
#' used by `deliverables/RDT-20/logs/oracle-recovery.json`.
#'
#' Requires `jsonlite`; declared in `Suggests`.
#'
#' @param panel A simulator return list (with `oracle_panel_att`).
#' @param synth_control_fit A return list from [synth_control_baci()].
#' @param paired_did_fit A return list from [baci_paired_did()] (or
#'   `NULL` to skip the DiD block).
#' @param cluster_boot A return list from
#'   [proponent_cluster_bootstrap()] (or `NULL`).
#' @param path Character. Output JSON path.
#' @return Invisibly returns `path`. Side effect: writes the log.
#' @export
write_oracle_recovery_log <- function(panel,
                                      synth_control_fit,
                                      paired_did_fit = NULL,
                                      cluster_boot = NULL,
                                      path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("write_oracle_recovery_log: package 'jsonlite' is required.")
  }
  stopifnot(
    is.list(panel), "oracle_panel_att" %in% names(panel),
    is.list(synth_control_fit),
    "panel_att_hat" %in% names(synth_control_fit),
    is.character(path), length(path) == 1L
  )

  log_entry <- list(
    generated_at      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    oracle_panel_att  = panel$oracle_panel_att,
    sc_panel_att_hat  = synth_control_fit$panel_att_hat,
    sc_minus_oracle   = synth_control_fit$panel_att_hat -
      panel$oracle_panel_att,
    sc_pre_fit_rmse   = synth_control_fit$pre_fit_rmse
  )
  if (!is.null(paired_did_fit)) {
    log_entry$did_att_hat <- paired_did_fit$att_hat
    log_entry$did_se_cluster <- paired_did_fit$se_cluster
    log_entry$did_minus_oracle <- paired_did_fit$att_hat -
      panel$oracle_panel_att
  }
  if (!is.null(cluster_boot)) {
    log_entry$bootstrap_estimator <- cluster_boot$estimator
    log_entry$bootstrap_se <- cluster_boot$se
    log_entry$bootstrap_ci_lo <- cluster_boot$ci_lo
    log_entry$bootstrap_ci_hi <- cluster_boot$ci_hi
    log_entry$bootstrap_n_boot <- cluster_boot$n_boot
  }
  log_entry$session_info <- utils::capture.output(
    print(utils::sessionInfo())
  )

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(jsonlite::toJSON(log_entry, dataframe = "rows",
                              auto_unbox = TRUE, pretty = TRUE,
                              null = "null", digits = NA), con = con)
  invisible(path)
}
