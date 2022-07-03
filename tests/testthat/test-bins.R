
test_that("Binning works", {
  expect_error(get_bins(exposure = 1:100, outcome = 100:1), NA)
  expect_error(get_bins(list(exposure = 1:100, outcome = 100:1)), NA)

  test_bins <- get_bins(exposure = 1:100, outcome = 100:1)
  expect_equal(test_bins[, -c("exposure_bin")],
               suppressWarnings(
                 data.table(exposure_median = 1:100,
                            outcome_median = 100:1,
                            outcome_mean = 100:1,
                            outcome_sd = numeric())
               ))
  expect_equal(get_bins(exposure = 1:100, outcome = 100:1, bins = 1e5),
               test_bins)

  test_bins <- get_bins(exposure = 1:100, outcome = 100:1, bins = 10)
  expect_equal(test_bins[, -c("exposure_bin")],
               suppressWarnings(
                 data.table(exposure_median = seq(5.5, 95.5, by = 10),
                            outcome_median = seq(95.5, 5.5, by = -10),
                            outcome_mean = seq(95.5, 5.5, by = -10),
                            outcome_sd = rep(sd(1:10), 10))
               ))
})
