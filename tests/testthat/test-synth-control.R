test_that("project_simplex returns a valid simplex projection", {
  set.seed(1L)
  v <- stats::rnorm(8)
  p <- bacipair:::project_simplex(v)
  expect_equal(sum(p), 1, tolerance = 1e-10)
  expect_true(all(p >= 0))
  ## Already on simplex -> identity (within tolerance).
  q <- rep(1 / 5, 5)
  expect_equal(bacipair:::project_simplex(q), q, tolerance = 1e-10)
})

test_that("solve_simplex_lsq recovers the noiseless solution", {
  set.seed(2L)
  X <- matrix(stats::rnorm(50 * 6), nrow = 50, ncol = 6)
  w_true <- bacipair:::project_simplex(stats::runif(6))
  y <- as.numeric(X %*% w_true)
  fit <- bacipair:::solve_simplex_lsq(X, y, max_iter = 5000L,
                                      tol = 1e-14)
  expect_equal(sum(fit$weights), 1, tolerance = 1e-8)
  expect_true(fit$objective < 1e-6)
})

test_that("synth_control_baci recovers the oracle within kill-gate floor", {
  sim <- cached_sim(n_pairs = 60L)
  fit <- bacipair::synth_control_baci(sim$panel)
  expect_true(is.numeric(fit$panel_att_hat))
  ## Operational floor (RDT-8 §5): aggregate bias <= 0.5 pp on a single
  ## small-n simulation. The tighter 0.15 pp gate uses ~12 replicates.
  expect_lt(abs(fit$panel_att_hat - sim$oracle_panel_att), 0.5)
})

test_that("synth_control_baci validates inputs", {
  sim <- cached_sim()
  expect_error(
    bacipair::synth_control_baci(as.data.frame(sim$panel)),
    regexp = "is.data.table"
  )
  small_panel <- sim$panel[get("pair_id") <= 4]
  expect_error(
    bacipair::synth_control_baci(small_panel),
    regexp = "donor pool"
  )
})

test_that("synth_control_baci runs on a real panel without oracle columns (F1)", {
  sim <- cached_sim(n_pairs = 60L)
  real_panel <- data.table::copy(sim$panel)
  real_panel[, c("true_additionality_max", "has_real_additionality",
                 "ite") := NULL]
  fit <- bacipair::synth_control_baci(real_panel)
  expect_true(is.numeric(fit$panel_att_hat))
  expect_false("oracle_unit" %in% names(fit$per_unit_att))
})

test_that("synth_control_baci abstains with a classed condition on a small donor pool (F4)", {
  sim <- cached_sim()
  small_panel <- sim$panel[get("pair_id") <= 4]
  cond <- tryCatch(
    bacipair::synth_control_baci(small_panel),
    error = function(e) e
  )
  expect_s3_class(cond, "bacipair_abstention")
  expect_true(any(grepl("_refusal$", class(cond))))
  expect_match(conditionMessage(cond), "donor pool")
})

test_that("synth_control_baci abstains with a classed condition on too few pre-years (F4)", {
  sim <- cached_sim(n_pairs = 30L)
  ## The default simulator panel's earliest-treated pairs (treatment_year
  ## 2013) have exactly 3 pre-treatment donor years (2010-2012); raising
  ## min_pre_years to 4 puts them below the floor without needing to
  ## mutate the panel.
  cond <- tryCatch(
    bacipair::synth_control_baci(sim$panel, min_pre_years = 4L),
    error = function(e) e
  )
  expect_s3_class(cond, "bacipair_abstention")
  expect_true(any(grepl("_refusal$", class(cond))))
})
