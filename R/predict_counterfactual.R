#' Counterfactual Prediction
#' @description Obtain counterfactual prediction of a fit.
#'
#' @param fit.j fitted object for trt j.
#' @param fit.k fitted object for trt k.
#' @param treatment name of treatment column
#' @param treatments_for_compare (`character`) Treatments for comparison
#' @param data (`data.frame`) raw dataset.
#' @param prob_mat (`data.frame`) treatment assignment probabilities
#' @param stabilize stabilize
#'
#' @return Numeric matrix of counter factual prediction.
#'
#' @export
predict_counterfactual <- function(fit.j,fit.k, treatment, treatments_for_compare, prob_mat, data,stabilize) {
  UseMethod("predict_counterfactual", fit.j)
}

#' @export
predict_counterfactual.lm <- function(fit.j,fit.k, treatment, treatments_for_compare,
                                      prob_mat, data = merge(find_data(fit.j),find_data(fit.k)),stabilize) {
  checkmate::assert_data_frame(data)
  checkmate::assert_subset(treatment, colnames(data))
  formula <- formula(fit.j)

  treatments_for_compare <- as.factor(treatments_for_compare)
  data[[treatment]] <- as.factor(data[[treatment]])
  checkmate::assert(
    checkmate::test_character(data[[treatment]]),
    checkmate::test_factor(data[[treatment]])
  )

  ECE_size = nrow(data)

  trt_lvls <- treatments_for_compare
  n_lvls <- length(treatments_for_compare)

  preds_linear.j <- predict(fit.j, newdata = data, type="link")
  preds_linear.k <- predict(fit.k, newdata = data, type="link")

  preds.j <- predict(fit.j, newdata = data, type="response")
  preds.k <- predict(fit.k, newdata = data, type="response")

  pred_linear <- matrix(c(preds_linear.j, preds_linear.k), ncol = n_lvls, dimnames = list(row.names(data), trt_lvls))
  ret <- matrix(c(preds.j, preds.k), ncol = n_lvls, dimnames = list(row.names(data), trt_lvls))
  y <- data[[all.vars(fit.j$formula)[1]]]

  estimation <- estimation_wt(ret, y, treatment, treatments_for_compare, data, prob_mat, stabilize=stabilize)

  structure(
    .Data = estimation,
    sample_size = ECE_size,
    predictions = ret,
    predictions_linear = pred_linear,
    response = y,
    fit.j = fit.j,
    fit.k = fit.k,
    prob_mat = prob_mat[treatments_for_compare],
    treatment_name = treatment,
    treatment = data[[treatment]],
    class = "prediction_cf"
  )
}

#' @export
predict_counterfactual.glm <- function(fit.j,fit.k, treatment, treatments_for_compare,
                                       prob_mat, data = merge(find_data(fit.j),find_data(fit.k)),stabilize, ...) {
  predict_counterfactual.lm(fit.j = fit.j,fit.k=fit.k, treatment = treatment,
                            treatments_for_compare = treatments_for_compare,prob_mat, data = data,stabilize,...)
}
