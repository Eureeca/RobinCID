#' S3 Methods for `prediction_cf`
#' @param x (`prediction_cf`)\cr the obtained counter-factual prediction object.
#' @name prediction_cf_methods
#' @return No return value, called for side effects
NULL

#' @describeIn prediction_cf_methods prints the prediction_cf object.
#' @exportS3Method
#' @keywords internal
print.prediction_cf <- function(x, ...) {
  cat("counter-factual prediction\n")
  print(x)
  cat("\n")
}
