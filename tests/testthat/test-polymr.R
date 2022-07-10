set.seed(1234)
my_data <- new_PolyMRDataSim()
polymr_model <- new_PolyMRModel(my_data)
eo_model <- new_EOModel(my_data)
test_that("Increasing the order of the function works", {
  expect_equal(increase_order_of_outcome_model(polymr_model),
               update_outcome_model(polymr_model, 1:3))
  expect_equal(increase_order_of_outcome_model(eo_model),
               update_outcome_model(eo_model, 1:3))

  expect_equal(increase_order_of_outcome_model(polymr_model, power_step = 3),
               update_outcome_model(polymr_model, 1:4, 1:4))
  expect_equal(increase_order_of_outcome_model(eo_model, power_step = 3),
               update_outcome_model(eo_model, 1:4))

  expect_equal(update_outcome_model(polymr_model, 1:2) |>
                 increase_order_of_outcome_model(),
               update_outcome_model(polymr_model, 1:4))
  expect_equal(update_outcome_model(eo_model, 1:2) |>
                 increase_order_of_outcome_model(),
               update_outcome_model(eo_model, 1:4))
})

test_that("Dropping the weakest term works", {
  expect_equal(polymr_model |>
                 update_outcome_model(1:3) |>
                 drop_least_significant_term(),
               polymr_model |>
                 update_outcome_model(1:2))
  expect_equal(eo_model |>
                 update_outcome_model(1:3) |>
                 drop_least_significant_term(),
               eo_model |>
                 update_outcome_model(1:2))

  expect_equal(polymr_model |>
                 update_outcome_model(1:3) |>
                 drop_least_significant_term(drop_higher_control_function_powers = FALSE),
               polymr_model |>
                 update_outcome_model(1:2, 1:3))
})

test_that("Ending up with an empty model is handled properly", {
  null_data <- new_PolyMRDataSim(causal_function = \(x) 0)
  null_model <- new_PolyMRModel(null_data)
  expect_warning(null_model |>
                   select_best_outcome_model(p_thr_drop = 1e-3),
                 "No significant terms")
  expect_equal(suppressWarnings(
                  select_best_outcome_model(null_model, p_thr_drop = 1e-3)
               )$outcome_model |>
                 coef(),
               get_null_model(null_model) |>
                 coef())
  expect_warning(null_model |>
                   update_outcome_model(1:4) |>
                   select_best_outcome_model(p_thr_drop = 1e-3),
                 "No significant terms")
})

test_that("Selecting the best outcome model works", {
  expect_equal(select_best_outcome_model(polymr_model,
                                         p_thr_add = .025,
                                         p_thr_drop = 1e-3),
               update_outcome_model(polymr_model, 1:2))
  expect_equal(select_best_outcome_model(eo_model,
                                         p_thr_add = .025,
                                         p_thr_drop = 1e-3),
               update_outcome_model(eo_model, 1:2))

  expect_equal(update_outcome_model(polymr_model, 1:10) |>
                 select_best_outcome_model(p_thr_drop = 1e-3),
               update_outcome_model(polymr_model, 1:2))
  expect_equal(update_outcome_model(eo_model, 1:10) |>
                 select_best_outcome_model(p_thr_drop = 1e-3),
               update_outcome_model(eo_model, 1:2))
})

test_that("polymr() returns the appropriate warnings", {
  expect_warning(polymr(my_data$exposure,
                        my_data$outcome,
                        my_data$genotypes,
                        max_control_function_power = 2))

  full_polymr <- polymr(my_data$exposure,
                        my_data$outcome,
                        my_data$genotypes,
                        p_thr_add = .025,
                        p_thr_drop = 1e-3)
  expect_equal(full_polymr$polymr,
               select_best_outcome_model(polymr_model |> update_outcome_model(1:10),
                                         p_thr_add = .025,
                                         p_thr_drop = 1e-3) |>
                 return_model())
})













