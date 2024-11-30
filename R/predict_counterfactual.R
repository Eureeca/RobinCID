#' Counterfactual Prediction
#' @description Obtain counterfactual prediction of a fit.
#'
#' @param fit fitted object.
#' @param treatment name of treatment column
#' @param treatments_for_compare (`character`) Treatments for comparison
#' @param data (`data.frame`) raw dataset.
#'
#' @return Numeric matrix of counter factual prediction.
#'
#' @export
predict_counterfactual <- function(fit, treatment, treatments_for_compare, prob_mat, data,...) {
  UseMethod("predict_counterfactual", fit)
}

#' @export
predict_counterfactual.lm <- function(fit, treatment, treatments_for_compare, prob_mat, data = find_data(fit), ...) {
  assert_data_frame(data)
  assert_subset(treatment, colnames(data))
  formula <- formula(fit)
  assert_subset(treatment, all.vars(formula[[3]]))

  treatments_for_compare <- as.factor(treatments_for_compare)
  data[[treatment]] <- as.factor(data[[treatment]])
  assert(
    test_character(data[[treatment]]),
    test_factor(data[[treatment]])
  )


  trt_lvls <- levels(treatments_for_compare)
  n_lvls <- length(trt_lvls)

  df <- lapply(
    data,
    function(i) {
      rep(i, times = n_lvls)
    }
  )

  df[[treatment]] <- rep(trt_lvls, each = nrow(data))

  mm <- model.matrix(fit, data = df)
  pred_linear <- mm %*% coefficients(fit)
  preds <- family(fit)$linkinv(pred_linear)

  ret <- matrix(preds, ncol = n_lvls, dimnames = list(row.names(data), trt_lvls))
  y <- model.response(fit$model)
  residual <- y - fitted(fit)

  estimation <- estimation_wt(ret, y, treatment, treatments_for_compare, data, prob_mat, stabilize=stabilize)

  structure(
    .Data = estimation,
    residual = residual,
    predictions = ret,
    predictions_linear = pred_linear,
    response = y,
    fit = fit,
    model_matrix = mm,
    treatment_name = treatment,
    treatment = data[[treatment]],
    class = "prediction_cf"
  )
}

#' @export
predict_counterfactual.glm <- function(fit, treatment, treatments_for_compare, prob_mat, data = find_data(fit), ...) {
  predict_counterfactual.lm(fit = fit, treatment = treatment,
                            treatments_for_compare = treatments_for_compare,prob_mat, data = data,...)
}
