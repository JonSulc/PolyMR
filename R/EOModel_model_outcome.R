model_outcome <- function(eo_model){
  eo_model <- create_outcome_predictors_table(eo_model)
  eo_model$outcome_model <-
    lm(eo_model$outcome ~ ., data = eo_model$outcome_predictors)
  eo_model <- calculate_summary_statistics(eo_model)
  eo_model
}


create_outcome_predictors_table <- function(x, ...){
  UseMethod("create_outcome_predictors_table")
}
create_outcome_predictors_table.EOModel <- function(eo_model){
  if (length(eo_model$exposure_powers) > 0){
    exposure_poly <-
      outer(eo_model$exposure,
            eo_model$exposure_powers,
            "^")
    exposure_poly <- data.table::as.data.table(exposure_poly)
    colnames(exposure_poly) <- paste0("x", eo_model$exposure_powers)
  } else
    exposure_poly <- NULL

  if (length(eo_model$control_function_powers) > 0){
    control_function_poly <-
      outer(eo_model$control_function,
            eo_model$control_function_powers,
            "^")
    control_function_poly <- data.table::as.data.table(control_function_poly)
    colnames(control_function_poly) <- paste0("ex", eo_model$control_function_powers)
  } else
    control_function_poly <- NULL

  eo_model$outcome_predictors <- cbind(exposure_poly,
                                       control_function_poly)
  eo_model
}


update_outcome_model <- function(x, ...){
  UseMethod("update_outcome_model")
}

update_outcome_model.EOModel <- function(eo_model,
                                         exposure_powers){
  eo_model$exposure_powers <- exposure_powers
  model_outcome(eo_model)
}
