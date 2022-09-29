set.seed(1234)

test_that("Plotting works in conjunction with term filtering", {
  my_data <- PolyMR:::new_PolyMRDataSim(
    sample_size = 10000,
    causal_function = function(x) 0.05 * exp(x)
  )

  expect_error(
    with(my_data,
         polymr(exposure, outcome, genotypes)) |>
      plot_polymr(),
    NA
  )

  expect_error(
    with(my_data,
         polymr(exposure, outcome, genotypes,
                p_thr_drop = NULL)) |>
      plot_polymr(),
    NA
  )
  predict.EOModel
})
