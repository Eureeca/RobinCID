#' Inverse probability weighting based model
#'
#' @param formula (`formula`) A formula of analysis.
#' @param data (`data.frame`) Input data frame.
#' @param treatment (`character(1)`) A string name of treatment assignment.
#' @param prob_mat (`data.frame`) Treatment assignment probability matrix
#' @param treatments_for_compare (`vector`) Treatments for comparison
#' @param contrast (`function` or `character(1)`) A function to calculate the treatment effect, or character of
#' "difference", "risk_ratio", "odds_ratio" for default contrasts.
#' @param contrast_jac (`function`) A function to calculate the Jacobian of the contrast function. Ignored if using
#' default contrasts.
#' @param family (`family`) A family object of the glm model.
#' @param stabilize (`logical`) Whether to stabilize
#' @param alpha (`double`) Nominal level
#' @param method Estimation method. Either "wt" or "ps".
#' @param ... Additional arguments passed to `glm`
#'
#' @export
#'
#' @examples
#' robin_wt(
#'   formula = y ~ xb + xc,
#'   data = example,
#'   treatment = "treatment",
#'   prob_mat = prob_mat,
#'   treatments_for_compare = c("1","2"),
#'   contrast = "difference",
#'   contrast_jac=NULL,
#'   family=gaussian(),
#'   stabilize=TRUE,
#'   alpha=0.05)
robin_estimate <- function(formula, data, treatment, prob_mat, treatments_for_compare,
                     contrast = "difference", contrast_jac=NULL, family=gaussian(), stabilize=TRUE, alpha=0.05,method=NULL,...){

  # Example: robin_wt(formula=y ~ treatment * x, data = dummy_data, treatment = treatment, prob_mat,
  # treatments_for_compare=, contrast = "difference", contrast_jac=NULL, family=gaussian(), stabilize=T)

  # Output:
  # Model        :  y ~ treatment * s1
  # Family: gaussian
  # Randomization probabilities (among the entire concurrent and eligible (ECE) samples):
  # unique levels of (pi_j(Z), pi_k(Z)) and sample sizes. Also include proportions and total size.
  #
  # Marginal Mean:
  #   Estimate Std.Err lower.CL upper.CL
  # trt_j         0.541      0.103    XX           XX
  # trt_k          0.725     0.103     XX           XX
  #
  # Treatment Effect:
  #   Contrast:  difference
  # Estimate Std.Err z.value lower.CL upper.CL Pr(>|z|)
  # trt_j - trt_k     0.541      0.103   5.26     XX           XX          X
  # ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

  checkmate::assert_subset(all.vars(formula), names(data))
  checkmate::assert_subset(treatment, names(data))
  checkmate::assert_subset(treatments_for_compare, names(prob_mat))
  checkmate::assert_subset(as.character(treatments_for_compare), as.character(data[[treatment]]))

  data[[treatment]] <- as.factor(data[[treatment]])
  treatments_for_compare <- factor(treatments_for_compare, levels=levels(data[[treatment]]))
  checkmate::assert(
    checkmate::test_character(data[[treatment]]),
    checkmate::test_factor(data[[treatment]])
  )

  # ECE sample
  data <- data[apply(prob_mat[treatments_for_compare]>0, 1, all),]
  prob_mat <- prob_mat[apply(prob_mat[treatments_for_compare]>0, 1, all),]

  fit.j <-  glm(formula, family = family, data = data[data[[treatment]]==treatments_for_compare[1],], ...)
  fit.k <-  glm(formula, family = family, data = data[data[[treatment]]==treatments_for_compare[2],], ...)
  pc <- predict_counterfactual(fit.j = fit.j, fit.k = fit.k, treatment = treatment,
                               treatments_for_compare = treatments_for_compare,
                               prob_mat = prob_mat, data = data,stabilize = stabilize, method = method)

  # has_interaction <- h_interaction(formula, treatment)


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
#' A wrapper function for `robin_estimate` with the method set to "wt".
#'
#' @inheritParams robin_estimate
#' @return See `robin_estimate` for details.
#' @export
#' @examples
#' robin_wt(
#'   formula = y ~ xb + xc,
#'   data = example,
#'   treatment = "treatment",
#'   prob_mat = prob_mat,
#'   treatments_for_compare = c("1", "2"),
#'   contrast = "difference",
#'   contrast_jac = NULL,
#'   family = gaussian(),
#'   stabilize = TRUE,
#'   alpha = 0.05
#' )
robin_wt <- function(...) {
  robin_estimate(method = "wt", ...)
}

#' Post-Stratification Based Inference
#'
#' A wrapper function for `robin_estimate` with the method set to "ps".
#'
#' @inheritParams robin_estimate
#' @return See `robin_estimate` for details.
#' @export
#' @examples
#' robin_ps(
#'   formula = y ~ xb + xc,
#'   data = example,
#'   treatment = "treatment",
#'   prob_mat = prob_mat,
#'   treatments_for_compare = c("1", "2"),
#'   contrast = "difference",
#'   contrast_jac = NULL,
#'   family = gaussian(),
#'   stabilize = TRUE,
#'   alpha = 0.05
#' )
robin_ps <- function(...) {
  robin_estimate(method = "ps", ...)
}
