my_data <- new_PolyMRDataSim()
test_that("PolyMRModel instanciation works", {
  expect_error(new_PolyMRModel(exposure = my_data$exposure,
                               outcome = my_data$outcome,
                               genotypes = my_data$genotypes),
               NA)
  expect_error(new_PolyMRModel(my_data), NA)

  expect_error(new_PolyMRModel(), "must all be provided")
  expect_error(new_PolyMRModel(exposure = my_data$exposure,
                               outcome = my_data$outcome),
               "must all be provided")
  expect_error(new_PolyMRModel(exposure = my_data$exposure[1:100],
                               outcome = my_data$outcome,
                               genotypes = my_data$genotypes),
               "sample sizes must match")
  expect_error(new_PolyMRModel(exposure = my_data$exposure,
                               outcome = my_data$outcome[1:100],
                               genotypes = my_data$genotypes),
               "sample sizes must match")
  expect_error(new_PolyMRModel(exposure = my_data$exposure,
                               outcome = my_data$outcome,
                               genotypes = my_data$genotypes[1:100, ]),
               "sample sizes must match")
  expect_error(new_PolyMRModel(exposure = my_data$exposure,
                               outcome = my_data$outcome,
                               genotypes = my_data$genotypes[ , 1:10]),
               NA)

  expect_equal(new_PolyMRModel(exposure = my_data$exposure,
                               outcome = my_data$outcome,
                               genotypes = my_data$genotypes),
               new_PolyMRModel(my_data))
})

test_that("EOModel instanciation works", {
  expect_error(new_EOModel(exposure = my_data$exposure,
                           outcome = my_data$outcome),
               NA)
  expect_error(new_EOModel(my_data), NA)

  expect_error(new_EOModel(), "[eE]xposure.+outcome")
  expect_error(new_EOModel(exposure = my_data$exposure),
               "[eE]xposure.+outcome")
  expect_error(new_EOModel(exposure = my_data$exposure[1:100],
                           outcome = my_data$outcome),
               "sample sizes must match")
  expect_error(new_EOModel(exposure = my_data$exposure,
                           outcome = my_data$outcome[1:100]),
               "sample sizes must match")

  expect_equal(new_EOModel(exposure = my_data$exposure,
                           outcome = my_data$outcome),
               new_EOModel(my_data))
})

polymr_model <- new_PolyMRModel(my_data)
test_that("PolyMRModel genotype-exposure association works", {
  expect_equal(polymr_model$exposure_model$coefficients |>
                 unname(),
               lm(scale(my_data$exposure) ~ my_data$genotypes)$coefficients |>
                 unname())
})

test_that("PolyMRModel control function calculation works", {
  expect_equal(c(polymr_model$genotypes %*% polymr_model$exposure_model$coef[-1]
                 + residuals(polymr_model$exposure_model)
                 + polymr_model$exposure_model$coef[1]),
               polymr_model$exposure)
})

eo_model <- new_EOModel(my_data)
test_that("Regression table works", {
  expect_identical(create_outcome_predictors_table(polymr_model)$outcome_predictors,
                   data.table::data.table(x1 = polymr_model$exposure,
                                          ex1 = residuals(polymr_model$exposure_model)))
  expect_identical(create_outcome_predictors_table(eo_model)$outcome_predictors,
                   data.table::data.table(x1 = eo_model$exposure))

  expect_identical(update_outcome_model(polymr_model, 1:2)$outcome_predictors,
                   data.table::data.table(x1 = polymr_model$exposure,
                                          x2 = polymr_model$exposure^2,
                                          ex1 = residuals(polymr_model$exposure_model),
                                          ex2 = residuals(polymr_model$exposure_model)^2))
  expect_identical(update_outcome_model(eo_model, 1:2)$outcome_predictors,
                   data.table::data.table(x1 = eo_model$exposure,
                                          x2 = eo_model$exposure^2))

  expect_identical(update_outcome_model(polymr_model, 2)$outcome_predictors,
                   data.table::data.table(x2 = polymr_model$exposure^2,
                                          ex1 = residuals(polymr_model$exposure_model),
                                          ex2 = residuals(polymr_model$exposure_model)^2))
  expect_identical(update_outcome_model(eo_model, 2)$outcome_predictors,
                   data.table::data.table(x2 = eo_model$exposure^2))

  expect_identical(update_outcome_model(polymr_model, 2, 3:4)$outcome_predictors,
                   data.table::data.table(x2 = polymr_model$exposure^2,
                                          ex3 = residuals(polymr_model$exposure_model)^3,
                                          ex4 = residuals(polymr_model$exposure_model)^4))
})

test_that("Modeling the outcome works", {
  expect_error(model_outcome(polymr_model), NA)
  expect_error(model_outcome(eo_model), NA)

  expect_identical(
    coef(polymr_model$outcome_model),
    coef(lm(polymr_model$outcome ~ ., data = polymr_model$outcome_predictors))
  )
  expect_identical(
    coef(eo_model$outcome_model),
    coef(lm(eo_model$outcome ~ ., data = eo_model$outcome_predictors))
  )

  expect_warning(update_outcome_model(polymr_model, 12))

  expect_identical(
    coef(update_outcome_model(polymr_model, exposure_powers = 2, control_function_powers = 1:2)$outcome_model),
    coef(lm(polymr_model$outcome ~ .,
            data = data.table::data.table(x2 = polymr_model$exposure^2,
                                          ex1 = residuals(polymr_model$exposure_model),
                                          ex2 = residuals(polymr_model$exposure_model)^2)))
  )
  expect_identical(
    coef(update_outcome_model(eo_model, exposure_powers = 2)$outcome_model),
    coef(lm(eo_model$outcome ~ .,
            data = data.table::data.table(x2 = eo_model$exposure^2)))
  )
})

test_that("Auto-updating control function powers works", {
  expect_equal(update_control_function_powers(polymr_model),
               polymr_model)
  expect_equal(polymr_model |>
                 update_outcome_model(1:3),
               polymr_model |>
                 update_outcome_model(1:3, 1:3))
  expect_equal(polymr_model |>
                 update_outcome_model(1:3),
               polymr_model |>
                 update_outcome_model(1:3, 3) |>
                 update_outcome_model(1:3))
})

test_that("Creating the null model works", {
  expect_identical(
    coef(get_null_model(polymr_model)),
    coef(lm(polymr_model$outcome ~ ex1, data = polymr_model$outcome_predictors))
  )
  expect_identical(
    update_outcome_model(polymr_model, exposure_powers = 1:3, control_function_powers = 4) |>
      get_null_model() |>
      coef(),
    coef(lm(polymr_model$outcome ~ ex4,
            data = data.table::data.table(ex4 = polymr_model$outcome_predictors$ex1^4)))
  )
})

test_that("Creating the linear model works", {
  expect_identical(
    coef(get_linear_model(polymr_model)),
    coef(lm(polymr_model$outcome ~ x1 + ex1, data = polymr_model$outcome_predictors))
  )
  expect_identical(
    coef(get_linear_model(eo_model)),
    coef(eo_model$outcome_model)
  )

  expect_identical(
    update_outcome_model(polymr_model, exposure_powers = c(1,3), control_function_powers = 4) |>
      get_linear_model() |>
      coef(),
    coef(lm(polymr_model$outcome ~ x1 + ex4,
            data = data.table::data.table(
              x1 = polymr_model$outcome_predictors$x1,
              ex4 = polymr_model$outcome_predictors$ex1^4
            )))
  )
  expect_identical(
    update_outcome_model(eo_model, exposure_powers = c(1,3)) |>
      get_linear_model() |>
      coef(),
    coef(eo_model$outcome_model)
  )
})

test_that("Null model p-values work correctly", {
  expect_equal(
    get_pval_vs_null_model(polymr_model),
    0,
    tolerance = 1e-5
  )
  expect_equal(
    get_pval_vs_null_model(eo_model),
    0,
    tolerance = 1e-5
  )
})

test_that("Linear model p-values work correctly", {
  expect_equal(get_pval_vs_linear_model(polymr_model), 1)
  expect_equal(get_pval_vs_linear_model(eo_model), 1)

  quadratic_model <- update_outcome_model(polymr_model, 1:2)
  expect_equal(
    get_pval_vs_linear_model(quadratic_model),
    get_lrt_pval_alternate_model(quadratic_model,
                                 get_linear_model(quadratic_model))
  )
  quadratic_model <- update_outcome_model(eo_model, 1:2)
  expect_equal(
    get_pval_vs_linear_model(quadratic_model),
    get_lrt_pval_alternate_model(quadratic_model,
                                 get_linear_model(quadratic_model))
  )
})

test_that("Getting the p-values of the last-added terms works", {
  new_model <- increase_order_of_outcome_model(polymr_model)
  expect_equal(are_terms_significant(new_model, 3) |>
                 length(),
               3)
})

test_that("Model summary statistics are properly updated", {
  expect_equal(polymr_model$pval_linear_model, 1)
  expect_equal(eo_model$pval_linear_model, 1)

  expect_true(update_outcome_model(polymr_model, 1:2)$pval_linear_model < 1)
  expect_true(update_outcome_model(eo_model, 1:2)$pval_linear_model < 1)
})

test_that("Calculating the variance-covariance (vcov) matrix works", {
  expect_equal(calculate_vcov(eo_model)$vcov |> diag(),
               summary(eo_model$outcome_model)$coef[ , "Std. Error"]^2)
  expect_equal(calculate_vcov(polymr_model)$vcov |> diag(),
               summary(polymr_model$outcome_model)$coef[ , "Std. Error"]^2)
})

test_that("Cleanup works", {
  clean_eo_model <- cleanup(eo_model)
  expect_true(setequal(
    names(clean_eo_model),
    c("outcome_model",
      "pval_null_model",
      "pval_linear_model",
      "r_squared")
  ))
  expect_s3_class(clean_eo_model, c("EOModel"))
  expect_error(summary(clean_eo_model$outcome_model), NA)

  clean_polymr_model <- cleanup(polymr_model)
  expect_true(setequal(
    names(clean_polymr_model),
    c("outcome_model",
      "pval_null_model",
      "pval_linear_model",
      "r_squared")
  ))
  expect_s3_class(clean_polymr_model, c("PolyMRModel", "EOModel"))
  expect_error(summary(clean_polymr_model$outcome_model), NA)
})













