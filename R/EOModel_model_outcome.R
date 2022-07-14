model_outcome <- function(eo_model){
  eo_model <- create_outcome_predictors_table(eo_model)
  eo_model$outcome_model <-
    lm(eo_model$outcome ~ ., data = eo_model$outcome_predictors)
  eo_model <- eo_model |>
    calculate_summary_statistics()
  eo_model
}


create_outcome_predictors_table <- function(x, ...){
  UseMethod("create_outcome_predictors_table")
}
create_outcome_predictors_table.EOModel <- function(eo_model){
  if (length(eo_model$exposure_powers) > 0){
    eo_model$outcome_predictors <-
      outer(eo_model$exposure,
            eo_model$exposure_powers,
            "^") |>
      data.table::as.data.table()
    colnames(eo_model$outcome_predictors) <-
      paste0("x", eo_model$exposure_powers)
  } else
    eo_model$outcome_predictors <- NULL
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
