preprocess_data <- function(eo_model, ...){
  UseMethod("preprocess_data")
}
preprocess_data.EOModel <- function(eo_model){
  scale_phenotypes(eo_model)
}
scale_phenotypes <- function(eo_model){
  # scale(x)[,] used to drop attributes returned by scale()
  eo_model[c("exposure", "outcome")] <-
    lapply(eo_model[c("exposure", "outcome")], \(x) scale(x)[,])
  eo_model
}
