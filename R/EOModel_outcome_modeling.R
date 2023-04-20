#' @import data.table

model_outcome <- function(eo_model) {
  eo_model <- create_outcome_predictors_table(eo_model)
  eo_model$outcome_model <-
    lm(eo_model$outcome ~ ., data = eo_model$outcome_predictors)
  eo_model <- eo_model |>
    calculate_summary_statistics()
  eo_model
}


create_outcome_predictors_table <- function(x, ...) {
  UseMethod("create_outcome_predictors_table")
}
create_outcome_predictors_table.EOModel <- function(eo_model) {
  eo_model$outcome_predictors <-
    create_power_table(eo_model$exposure,
                       eo_model$exposure_powers,
                       "exposure")
  eo_model
}
create_power_table <- function(base, powers, name) {
  if (length(powers) == 0)
    return(NULL)
  power_table <- outer(base, powers, "^") |>
    data.table::as.data.table()
  colnames(power_table) <-
    paste0(name, powers)
  power_table
}


update_outcome_model <- function(x, ...){
  UseMethod("update_outcome_model")
}

update_outcome_model.EOModel <- function(eo_model,
                                         exposure_powers){
  eo_model$exposure_powers <- exposure_powers
  model_outcome(eo_model)
}
