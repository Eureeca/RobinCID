#' @noRd
robin_estimate <- function(formula, data, treatment, probabilities = NULL, post_strat = NULL, treatments_for_compare=NULL,
                     contrast = "difference", contrast_jac=NULL, family=gaussian(), stabilize=TRUE, alpha=0.05, method=NULL,...){

  assert_subset(all.vars(formula), names(data))
  assert_subset(treatment, names(data))
  assert_subset(as.character(treatments_for_compare), as.character(data[[treatment]]))

  if(identical(method, "ps")){
    check_result <- prob_strata_check(data, treatment, probabilities, post_strat, treatments_for_compare)
    prob_mat <- check_result$prob_mat
    post_strat <- check_result$post_strat
  } else {
    prob_mat <- data[probabilities]
  }

  data[[treatment]] <- as.factor(data[[treatment]])
  treatments_for_compare <- factor(treatments_for_compare, levels=levels(data[[treatment]]))
  assert(
    test_character(data[[treatment]]),
    test_factor(data[[treatment]])
  )

  # ECE sample
  ECE_subset <- apply(data[as.character(treatments_for_compare)] > 0, 1, all)
  data <- data[ECE_subset, ]
  if (nrow(data) < 1)
    stop("ECE sample size is 0!")
  prob_mat <- prob_mat[ECE_subset, ]

  if (identical(family$family, "Negative Binomial(NA)")) {
    fit.j <-  MASS::glm.nb(formula, family = family, data = data[data[[treatment]]==treatments_for_compare[1],], ...)
    fit.k <-  MASS::glm.nb(formula, family = family, data = data[data[[treatment]]==treatments_for_compare[2],], ...)
  } else {
    fit.j <-  glm(formula, family = family, data = data[data[[treatment]]==treatments_for_compare[1],], ...)
    fit.k <-  glm(formula, family = family, data = data[data[[treatment]]==treatments_for_compare[2],], ...)
  }
  pc <- predict_counterfactual(fit.j = fit.j, fit.k = fit.k, treatment = treatment,
                               treatments_for_compare = treatments_for_compare,
                               prob_mat = prob_mat, post_strat = post_strat,
                               data = data, stabilize = stabilize, method = method)


  if (identical(contrast, "difference")) {
    difference(pc, alpha = alpha)
  } else if (identical(contrast, "risk_ratio")) {
    risk_ratio(pc, alpha = alpha)
  } else if (identical(contrast, "odds_ratio")) {
    odds_ratio(pc, alpha = alpha)
  } else {
    assert_function(contrast)
    assert_function(contrast_jac, null.ok = TRUE)
    if (is.null(contrast_jac)) {
      contrast_jac <- function(x) {
        numDeriv::jacobian(contrast, x)
      }
    }
    treatment_effect(pc, eff_measure = contrast, eff_jacobian = contrast_jac, alpha = alpha)
  }
}

#' Inverse Probability Weighting Based Inference
#'
#' Provides robust inference methods via inverse probability weighting.
#'
#' @param formula (`formula`) A formula of analysis.
#' @param data (`data.frame`) Input data frame.
#' @param treatment (`character`) A string name of treatment assignment.
#' @param probabilities (`vector`) A character vector specifying column names of treatment assignment probabilities in the `data`.
#' The column names must represent the treatment levels in the study and should include at least the treatments specified in the `treatments_for_compare`.
#' @param treatments_for_compare (`vector`) A character vector specifying exactly two treatment levels to compare.
#' The specified treatments must be present in both the `treatment` variable of the `data` and the column names of the `probability`.
#' @param contrast (`function` or `character`) A function to calculate the treatment effect, or character of
#' "difference", "risk_ratio", "odds_ratio" for default contrasts.
#' @param contrast_jac (`function`) A function to calculate the Jacobian of the contrast function. Ignored if using
#' default contrasts.
#' @param family (`family`) A family object of the glm model. Default: `gaussian()`.
#' @param stabilize (`logical`) Whether to stabilize. Default: TRUE.
#' @param alpha (`double`) Nominal level. Default: 0.05.
#' @param ... Additional arguments passed to `glm`.
#' @details
#' If family is `MASS::negative.binomial(NA)`, the function will use `MASS::glm.nb` instead of `glm`.
#' @export
#' @examples
#' robin_wt(
#'   formula = y ~ xb + xc,
#'   data = example,
#'   treatment = "treatment",
#'   probabilities = c("trt.1", "trt.2", "trt.3", "trt.4"),
#'   treatments_for_compare = c("trt.1", "trt.2"),
#'   contrast = "difference",
#'   contrast_jac = NULL,
#'   family = gaussian(),
#'   stabilize = TRUE,
#'   alpha = 0.05
#' )
robin_wt <- function(formula, data, treatment, probabilities, treatments_for_compare,
                     contrast = "difference", contrast_jac=NULL, family=gaussian(), stabilize=TRUE, alpha=0.05,...) {
  if(is.null(probabilities)) stop("Assignment probabilities MUST be provided.")
  assert_subset(treatments_for_compare, probabilities)

  robin_estimate(formula = formula, data = data, treatment = treatment,
                 probabilities = probabilities, treatments_for_compare = treatments_for_compare,
                 contrast = contrast, contrast_jac = contrast_jac, stabilize = stabilize,
                 alpha = alpha, method = "wt", post_strat = NULL, ...)
}

#' Post-Stratification Based Inference
#'
#' Provides robust inference methods via post stratification.
#'
#' @inheritParams robin_wt
#' @param post_strat (`character`) A string name of post-stratification variable. Default: `NULL`
#' @details
#' If family is `MASS::negative.binomial(NA)`, the function will use `MASS::glm.nb` instead of `glm`.
#' @export
#' @examples
#' robin_ps(
#'   formula = y ~ xb + xc,
#'   data = example,
#'   treatment = "treatment",
#'   probabilities = c("trt.1", "trt.2", "trt.3", "trt.4"),
#'   post_strat = NULL,
#'   treatments_for_compare = c("trt.1", "trt.2"),
#'   contrast = "difference",
#'   contrast_jac = NULL,
#'   family = gaussian(),
#'   alpha = 0.05
#' )
robin_ps <- function(formula, data, treatment, probabilities = NULL, post_strat = NULL, treatments_for_compare,
                     contrast = "difference", contrast_jac = NULL, family = gaussian(), alpha = 0.05, ...) {

  robin_estimate(formula = formula, data = data, treatment = treatment,
                 probabilities = probabilities, post_strat = post_strat, treatments_for_compare = treatments_for_compare,
                 contrast = contrast, contrast_jac = contrast_jac, family=family,
                 alpha=alpha, method = "ps", ...)
}
