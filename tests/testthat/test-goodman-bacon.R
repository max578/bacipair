test_that("goodman_bacon_decomp returns a well-formed decomposition", {
  skip_if_not_installed("bacondecomp")
  sim <- cached_sim(n_pairs = 40L)
  gb <- bacipair::goodman_bacon_decomp(sim$panel)
  expect_named(gb, c("decomposition", "twfe_estimate",
                     "weight_summary"))
  expect_s3_class(gb$decomposition, "data.frame")
  expect_true(all(c("treated", "untreated", "estimate",
                    "weight", "type") %in% names(gb$decomposition)))
  expect_equal(sum(gb$decomposition$weight), 1, tolerance = 1e-6)
  expect_true(is.numeric(gb$twfe_estimate))
})

test_that("goodman_bacon_decomp errors on a malformed indicator", {
  skip_if_not_installed("bacondecomp")
  sim <- cached_sim()
  bad <- data.table::copy(sim$panel)
  bad[, "post_treatment" := as.numeric(get("post_treatment")) + 0.5]
  expect_error(
    bacipair::goodman_bacon_decomp(bad),
    regexp = "0/1 valued"
  )
})
