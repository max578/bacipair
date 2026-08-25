test_that("simulate_staggered_panel returns a well-formed object", {
  sim <- cached_sim()
  expect_named(sim, c("panel", "pair_table", "rainfall",
                      "oracle_panel_att", "oracle_unit_att",
                      "parameters"))
  expect_s3_class(sim$panel, "data.table")
  expect_true(all(c("pair_id", "proponent_id", "arm", "year",
                    "treatment_year", "post_treatment", "treated_arm",
                    "post_idx", "woody_cover", "ite") %in% names(sim$panel)))
  expect_equal(sort(unique(sim$panel$arm)), c("comparison", "credited"))
  expect_true(sim$oracle_panel_att > 0)
})

test_that("simulator is deterministic on the same seed", {
  sim_a <- bacipair::simulate_staggered_panel(seed = 20260520L,
                                              n_pairs = 20L,
                                              n_proponents = 6L)
  sim_b <- bacipair::simulate_staggered_panel(seed = 20260520L,
                                              n_pairs = 20L,
                                              n_proponents = 6L)
  expect_equal(sim_a$oracle_panel_att, sim_b$oracle_panel_att)
  expect_equal(sim_a$panel$woody_cover, sim_b$panel$woody_cover)
})

test_that("simulator validates input ranges", {
  expect_error(
    bacipair::simulate_staggered_panel(n_pairs = 1L),
    regexp = "n_pairs"
  )
  expect_error(
    bacipair::simulate_staggered_panel(t_min = 2025L, t_max = 2010L),
    regexp = "t_max"
  )
  expect_error(
    bacipair::simulate_staggered_panel(additionality_share = 1.5),
    regexp = "additionality_share"
  )
  expect_error(
    bacipair::simulate_staggered_panel(n_proponents = 0L,
                                       n_pairs = 10L),
    regexp = "n_proponents"
  )
})

test_that("oracle ATT is non-zero only on treated post-treatment rows", {
  sim <- cached_sim()
  zero_block <- sim$panel[treated_arm == 0L | post_treatment == 0L]
  expect_true(all(zero_block$ite == 0))
  treated_post <- sim$panel[treated_arm == 1L & post_treatment == 1L]
  ## With additionality_share = 0.18 some treated rows may have zero
  ## true_additionality_max; we only assert at least one positive ITE.
  expect_true(any(treated_post$ite > 0))
})
