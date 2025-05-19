#' Counterfactual Prediction
#' @description Obtain counterfactual prediction of a fit.
#'
#' @param fit.j fitted object for trt j.
#' @param fit.k fitted object for trt k.
#' @param treatment name of treatment column
#' @param treatments_for_compare (`character`) Treatments for comparison
#' @param data (`data.frame`) raw dataset.
#' @param prob_mat (`data.frame`) treatment assignment probabilities
#' @param post_strata (`character`) A string name of post-stratification variable.
#' @param stabilize stabilize
#' @param settings estimation setting
#'
#' @return A list of `prediction_cf` object with following elements:
#' - `estimation`: = estimation,
#' - `sample_size`: sample size of ECE population.
#' - `fit.j`: fitted model for the treatment specified in tx_to_compare.
#' - `fit.k`: fitted model for the treatment specified in tx_to_compare.
#'
#' @export
predict_counterfactual <- function(fit.j,fit.k, treatment, treatments_for_compare, prob_mat, post_strata, data, stabilize, settings) {
  UseMethod("predict_counterfactual", fit.j)
}

#' @export
predict_counterfactual.lm <- function(fit.j,fit.k, treatment, treatments_for_compare,
                                      prob_mat, post_strata, data = merge(find_data(fit.j),find_data(fit.k)),stabilize, settings) {
  checkmate::assert_data_frame(data)
  checkmate::assert_subset(treatment, colnames(data))
  formula <- formula(fit.j)

  treatments_for_compare <- as.character(treatments_for_compare)
  data[[treatment]] <- as.factor(data[[treatment]])
  checkmate::assert(
    checkmate::test_character(data[[treatment]]),
    checkmate::test_factor(data[[treatment]])
  )


  ECE_size <- nrow(data)

  trt_lvls <- treatments_for_compare
  n_lvls <- length(treatments_for_compare)

  preds_linear.j <- predict(fit.j, newdata = data, type="link")
  preds_linear.k <- predict(fit.k, newdata = data, type="link")

  preds.j <- predict(fit.j, newdata = data, type="response")
  preds.k <- predict(fit.k, newdata = data, type="response")

  pred_linear <- matrix(c(preds_linear.j, preds_linear.k), ncol = n_lvls, dimnames = list(row.names(data), trt_lvls))
  ret <- matrix(c(preds.j, preds.k), ncol = n_lvls, dimnames = list(row.names(data), trt_lvls))
  y <- data[[all.vars(fit.j$formula)[1]]]

  class(ret) <- settings$method
  estimation <- estimate_effect(ret, y, treatment, treatments_for_compare, data, prob_mat, post_strata, stabilize=stabilize)

  structure(
    list(
      estimation = estimation,
      sample_size = ECE_size,
      fit.j = fit.j,
      fit.k = fit.k
    ),
    predictions = ret,
    predictions_linear = pred_linear,
    response = y,
    prob_mat = prob_mat[treatments_for_compare],
    Z = attr(prob_mat, "Z"),
    post_strata = post_strata,
    treatment_name = treatment,
    treatment = data[[treatment]],
    data = data,
    class = "prediction_cf",
    settings = settings
  )
}

#' @export
predict_counterfactual.glm <- function(fit.j,fit.k, treatment, treatments_for_compare,
                                       prob_mat, post_strata, data = merge(find_data(fit.j),find_data(fit.k)),stabilize,settings){
  predict_counterfactual.lm(fit.j = fit.j,fit.k=fit.k, treatment = treatment,
                            treatments_for_compare = treatments_for_compare,prob_mat, post_strata = post_strata,
                            data = data, stabilize, settings)
}
