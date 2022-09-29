
my_data <- new_PolyMRDataSim(finalize = FALSE)
test_that("Simulated genotype generation works", {
  # skip("Time consuming")
  expect_error(new_PolyMRDataSim(finalize = FALSE), NA)
  expect_equal(length(my_data$exposure), 1e5)
  expect_equal(length(new_PolyMRDataSim(1e3, finalize = FALSE)$exposure),  1e3)
  expect_equal(dim(my_data$genotypes), c(1e5, 100))
  expect_equal(dim(new_PolyMRDataSim(sample_size = 1e4,
                                     n_exposure_snps = 35,
                                     finalize = FALSE)$genotypes),
               c(1e4, 35))

  expect_true(all(my_data$genotypes >= 0) && all(my_data$genotypes <= 2))
  expect_equal(apply(my_data$genotypes, 2, mean), 2*my_data$mafs, tolerance = .1)
  expect_equal(my_data$exposure_coefficients|>length(), my_data$n_exposure_snps)
  expect_equal(sum(my_data$exposure_coefficients^2), my_data$exposure_heritability)

  my_data <- new_PolyMRDataSim(sample_size = 10, n_exposure_snps = 1e4, finalize = FALSE)
  expect_true(all(colSums(my_data$genotypes) != 0))

  expect_equal(remove_invariant_snps(my_data), my_data)
})

test_that("Confounders work", {
  # skip("Time consuming")
  my_confounder <- new_Confounder(sample_size = 1e5)
  expect_equal(my_confounder$exposure_confounding_function(my_confounder$confounder_values),
               my_confounder$confounder_values * 0.2)
  expect_equal(my_confounder$outcome_confounding_function(my_confounder$confounder_values),
               my_confounder$confounder_values * 0.5)

  expect_equal(mean(my_confounder$confounder_values), 0)
  expect_equal(  sd(my_confounder$confounder_values), 1)
  expect_equal(length(my_confounder$confounder_values), 1e5)
})

my_data <- new_PolyMRDataSim(confounders_list = list(), finalize = FALSE)
test_that("Phenotype generation works", {
  # skip("Time consuming")
  exposure0 <- c(scale(my_data$genotypes)[,] %*% my_data$exposure_coefficients)
  expect_equal(my_data$exposure, exposure0)
  outcome0 <- my_data$causal_function(my_data$exposure)
  expect_equal(my_data$outcome, outcome0)

  my_confounder <- new_Confounder(1e5)
  my_data <- my_data + my_confounder
  exposure0 <- exposure0 + .2*my_confounder$confounder_values
  outcome0  <- my_data$causal_function(my_data$exposure) + .5*my_confounder$confounder_values
  expect_equal(length(my_data$confounders_list), 1)
  expect_equal(my_data$exposure, exposure0)
  expect_equal(my_data$outcome,  outcome0)

  my_confounder <- new_Confounder(
    sample_size = 1e5,
    exposure_confounding_function = get_polynomial_function(.1),
    outcome_confounding_function = get_polynomial_function(c(.05, .1))
  )
  my_data <- my_data + my_confounder
  exposure0 <- exposure0 + .1*my_confounder$confounder_values
  outcome0  <- my_data$causal_function(my_data$exposure) +
    .5*my_data$confounders_list[[1]]$confounder_values +
    .05*my_confounder$confounder_values +
    .1*my_confounder$confounder_values^2

  expect_equal(length(my_data$confounders_list), 2)
  expect_equal(my_data$exposure, exposure0)
  expect_equal(my_data$outcome,  outcome0)
})

test_that("Phenotype completion (adding error) works", {
  # skip("Time consuming")
  my_confounder <- new_Confounder(
    sample_size = 1e5,
    exposure_confounding_function = get_polynomial_function(.1),
    outcome_confounding_function = get_polynomial_function(c(.05, .1))
  )

  my_data <- my_data + my_confounder
  expect_true(var(my_data$exposure) < 1)
  expect_true(var(my_data$outcome)  < 1)
  my_data <- finalize_data(my_data)
  expect_equal(var(my_data$exposure), 1, tolerance = .01)
  expect_equal(var(my_data$outcome),  1, tolerance = .01)

  my_data <- new_PolyMRDataSim(finalize = TRUE)
  expect_equal(var(my_data$exposure), 1, tolerance = .01)
  expect_equal(var(my_data$outcome),  1, tolerance = .01)
})

test_that("Excess explained variance generates an error", {
  expect_error(new_PolyMRDataSim(causal_function = \(x) x),
               "variance of outcome exceeds 1")
  expect_error(new_PolyMRDataSim(exposure_heritability = .99),
               "variance.+exposure.+1")
})

test_that("GWS filtering SNPs works", {
  # skip("Time consuming")
  n_snps <- ncol(my_data$genotypes)
  my_data <- finalize_data(my_data)
  expect_true(ncol(my_data$genotypes) < n_snps)
  expect_error(validate_PolyMRDataSim(my_data), NA)
})

























