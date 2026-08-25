test_that("baci_paired_did returns a sensible point estimate and SE", {
  sim <- cached_sim(n_pairs = 60L)
  did <- bacipair::baci_paired_did(sim$panel)
  expect_type(did$att_hat, "double")
  expect_true(is.finite(did$att_hat))
  expect_true(did$se_cluster > 0)
  expect_equal(did$ci_hi - did$ci_lo,
               2 * stats::qnorm(0.975) * did$se_cluster,
               tolerance = 1e-10)
  expect_equal(did$n_pairs, 60L)
})

test_that("baci_paired_did is unbiased on a placebo (zero-effect) panel", {
  sim_placebo <- bacipair::simulate_staggered_panel(
    seed = 20260601L,
    n_pairs = 80L,
    additionality_share = 0.0
  )
  did <- bacipair::baci_paired_did(sim_placebo$panel)
  ## With zero true ATT and 80 pairs the CI should sit close to 0;
  ## allow generous tolerance to keep this stable across platforms.
  expect_lt(abs(did$att_hat), 0.5)
})

test_that("baci_paired_did errors on a missing arm level", {
  sim <- cached_sim()
  bad <- data.table::copy(sim$panel)
  bad <- bad[get("arm") == "comparison"]
  expect_error(
    bacipair::baci_paired_did(bad),
    regexp = "must contain"
  )
})

test_that("baci_paired_did errors on insufficient year coverage", {
  sim <- cached_sim()
  one_year <- data.table::copy(sim$panel[get("year") == 2015L])
  expect_error(
    bacipair::baci_paired_did(one_year),
    regexp = ">= 2 distinct years"
  )
})
