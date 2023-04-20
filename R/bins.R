#' @import data.table

get_bins <- function(model_data = NULL,
                     exposure = model_data$exposure,
                     outcome = model_data$outcome,
                     bins = 100){
  if ((is.null(model_data$exposure)  && is.null(exposure)) ||
      (is.null(model_data$outcome)   && is.null(outcome)))
    stop("Exposure and outcome must both be provided.")
  stopifnot(length(exposure) == length(outcome))

  if (bins > length(exposure))
    bins <- length(exposure)

  bin_boundaries <- quantile(exposure,
                             seq(0, 1, length.out = bins+1))

  eo_data <- data.table::data.table(
    exposure = exposure,
    outcome = outcome,
    exposure_bin = cut(exposure,
                       breaks = bin_boundaries,
                       include.lowest = TRUE,
                       ordered_result = TRUE)
  )

  eo_data[
    order(exposure_bin),
    .(exposure_median = median(exposure),
      outcome_median  = median(outcome),
      outcome_mean    = mean(outcome),
      outcome_sd      = sd(outcome)),
    by = exposure_bin
  ]
}



