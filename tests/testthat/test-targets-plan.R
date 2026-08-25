test_that("bacipair_targets_plan returns the expected node names", {
  skip_if_not_installed("targets")
  plan <- bacipair::bacipair_targets_plan()
  expect_named(plan, c("panel", "synth_control_fit", "paired_did_fit",
                       "cluster_boot", "goodman_bacon",
                       "oracle_recovery_log"))
  for (node in plan) {
    expect_s3_class(node, "tar_target")
  }
})

test_that("write_oracle_recovery_log writes a valid JSON log", {
  skip_if_not_installed("jsonlite")
  sim <- cached_sim(n_pairs = 25L)
  sc  <- bacipair::synth_control_baci(sim$panel)
  did <- bacipair::baci_paired_did(sim$panel)
  log_path <- file.path(tempdir(), "oracle-recovery-test.json")
  on.exit(unlink(log_path), add = TRUE)
  written <- bacipair::write_oracle_recovery_log(
    panel = sim,
    synth_control_fit = sc,
    paired_did_fit = did,
    cluster_boot = NULL,
    path = log_path
  )
  expect_equal(written, log_path)
  parsed <- jsonlite::fromJSON(log_path)
  expect_true(all(c("oracle_panel_att", "sc_panel_att_hat",
                    "did_att_hat") %in% names(parsed)))
  expect_equal(parsed$oracle_panel_att, sim$oracle_panel_att,
               tolerance = 1e-10)
})
