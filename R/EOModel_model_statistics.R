calculate_summary_statistics <- function(x, ...){
  UseMethod("calculate_summary_statistics")
}


calculate_summary_statistics.EOModel <- function(eo_model){
  eo_model$pval_null_model <- get_pval_vs_null_model(eo_model)
  eo_model$pval_linear_model <- get_pval_vs_linear_model(eo_model)
  eo_model$r_squared <- get_r_squared(eo_model)

  eo_model
}


.get_outcome_model_coefficients <- function(model,
                                            terms_to_include,
                                            values_to_return){
  all_coefficients <- model$outcome_model |>
    summary() |>
    coef()
  all_coefficients[terms_to_include, values_to_return]
}

get_outcome_model_coefficients <- function(x, ...){
  UseMethod("get_outcome_model_coefficients")
}

get_outcome_model_coefficients.EOModel <- function(
    eo_model,
    include_exposure = TRUE,
    include_intercept = FALSE,
    terms_to_include = NULL,
    values_to_return = c("Estimate", "Pr(>|t|)")
){
  if (include_exposure && length(eo_model$exposure_powers) > 0)
    terms_to_include <- c(paste0("x", eo_model$exposure_powers),
                          terms_to_include)
  if (include_intercept)
    terms_to_include <- c("(Intercept)", terms_to_include)

  .get_outcome_model_coefficients(eo_model,
                                  terms_to_include,
                                  values_to_return)
}

get_outcome_model_coefficients_exposure <- function(eo_model, ...){
  UseMethod("get_outcome_model_coefficients_exposure")
}

get_outcome_model_coefficients_exposure_pvalues <- function(eo_model){
  UseMethod("get_outcome_model_coefficients_exposure_pvalues")
}

get_outcome_model_coefficients_exposure_pvalues.EOModel <- function(eo_model){
  get_outcome_model_coefficients(eo_model,
                                 include_exposure = TRUE,
                                 values = "Pr(>|t|)")
}


get_linear_model <- function(eo_model){
  UseMethod("get_linear_model")
}

get_linear_model.EOModel <- function(eo_model){
  terms_to_keep <- colnames(eo_model$outcome_predictors)[
    grepl("^(x1$|[^x])", colnames(eo_model$outcome_predictors))
  ]
  outcome_predictors <- eo_model$outcome_predictors[ , ..terms_to_keep]

  if (!("x1" %in% colnames(eo_model$outcome_predictors)))
    outcome_predictors <- cbind(data.table::data.table(x1 = eo_model$exposure),
                                outcome_predictors)

  linear_model_terms <- colnames(outcome_predictors) |>
    paste(collapse = " + ")
  linear_model_formula <- paste(". ~", linear_model_terms) |>
    as.formula()

  update(eo_model$outcome_model, linear_model_formula, data = outcome_predictors)
}


get_lrt_pval_alternate_model <- function(eo_model, alternate_outcome_model){
  lrt <- lmtest::lrtest(eo_model$outcome_model, alternate_outcome_model)
  if (lrt[1, "LogLik"] < lrt[2, "LogLik"])
    return(1)
  lrt[2, "Pr(>Chisq)"]
}


get_pval_vs_null_model <- function(model, ...){
  UseMethod("get_pval_vs_null_model")
}

get_pval_vs_null_model.EOModel <- function(eo_model){
  model_fstatistics <- summary(eo_model$outcome_model)$fstatistic
  pf(q = model_fstatistics[['value']],
     df1 = model_fstatistics[['numdf']],
     df2 = model_fstatistics[['dendf']],
     lower.tail = FALSE)
}


get_pval_vs_linear_model <- function(model, ...){
  UseMethod("get_pval_vs_linear_model")
}

get_pval_vs_linear_model.EOModel <- function(eo_model){
  linear_model <- get_linear_model(eo_model)
  get_lrt_pval_alternate_model(eo_model, linear_model)
}


get_r_squared <- function(model){
  UseMethod("get_r_squared")
}

get_r_squared.EOModel <- function(eo_model){
  summary(eo_model$outcome_model)$r.squared
}


calculate_vcov <- function(eo_model){
  UseMethod("calculate_vcov")
}
calculate_vcov.EOModel <- function(eo_model){
  eo_model$vcov <- vcov(eo_model$outcome_model)
  eo_model
}
