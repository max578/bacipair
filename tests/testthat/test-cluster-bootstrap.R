test_that("proponent_cluster_bootstrap (paired_did) gives a sensible CI", {
  sim <- cached_sim(n_pairs = 40L)
  boot <- bacipair::proponent_cluster_bootstrap(
    sim$panel, estimator = "paired_did", n_boot = 50L,
    seed = 20260521L
  )
  expect_length(boot$draws, 50L)
  expect_true(boot$se > 0)
  expect_true(boot$ci_lo < boot$ci_hi)
  expect_equal(boot$estimator, "paired_did")
})

test_that("proponent_cluster_bootstrap is deterministic on the same seed", {
  sim <- cached_sim(n_pairs = 30L)
  a <- bacipair::proponent_cluster_bootstrap(
    sim$panel, estimator = "paired_did", n_boot = 50L, seed = 7L
  )
  b <- bacipair::proponent_cluster_bootstrap(
    sim$panel, estimator = "paired_did", n_boot = 50L, seed = 7L
  )
  expect_equal(a$draws, b$draws)
})

test_that("proponent_cluster_bootstrap rejects invalid n_boot", {
  sim <- cached_sim()
  expect_error(
    bacipair::proponent_cluster_bootstrap(
      sim$panel, estimator = "paired_did", n_boot = 5L
    ),
    regexp = "n_boot"
  )
})
