#' @noRd
robin_estimate <- function(formula, data, treatment,
                           probabilities = NULL, Z = NULL, Z_table = NULL,
                           post_strata = NULL, treatments_for_compare=NULL,
                           contrast = "difference", contrast_jac=NULL,
                           family=gaussian(), stabilize=TRUE, alpha=0.05, method=NULL, ...){

  assert_subset(all.vars(formula), names(data))
  assert_subset(treatment, names(data))
  assert_subset(as.character(treatments_for_compare), as.character(data[[treatment]]))


  if(!is.null(Z) & !is.null(Z_table)){
    data <- assign_prob_and_strata(data = data, Z = Z, Z_table = Z_table, probabilities = probabilities, merge_strata = T)
    post_strata <- "post_stratum"
  }

  if(identical(method, "ps")){
    check_result <- prob_strata_check(data, treatment, probabilities, post_strata, treatments_for_compare)
    prob_mat <- check_result$prob_mat
    post_strata <- check_result$post_strata
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
  attr(prob_mat, "Z") = Z

  if (identical(family$family, "Negative Binomial(NA)")) {
    fit.j <-  MASS::glm.nb(formula, family = family, data = data[data[[treatment]]==treatments_for_compare[1],], ...)
    fit.k <-  MASS::glm.nb(formula, family = family, data = data[data[[treatment]]==treatments_for_compare[2],], ...)
  } else {
    fit.j <-  glm(formula, family = family, data = data[data[[treatment]]==treatments_for_compare[1],], ...)
    fit.k <-  glm(formula, family = family, data = data[data[[treatment]]==treatments_for_compare[2],], ...)
  }
  pc <- predict_counterfactual(fit.j = fit.j, fit.k = fit.k, treatment = treatment,
                               treatments_for_compare = treatments_for_compare,
                               prob_mat = prob_mat, post_strata = post_strata,
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
#' @param Z (`vector`) The variable names in `data` and `Z_table` that define `stratum`.
#' @param Z_table (`data.frame`) A reference table where rows define strata and columns include `Z` and `probabilities`.
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
#' probabilities <- c("trt.1", "trt.2", "trt.3", "trt.4")
#'
#' robin_wt(
#'   formula = y ~ xb + xc,
#'   data = example,
#'   treatment = "treatment",
#'   probabilities = probabilities,
#'   treatments_for_compare = c("trt.1", "trt.2"),
#'   contrast = "difference",
#'   contrast_jac = NULL,
#'   family = gaussian(),
#'   stabilize = TRUE,
#'   alpha = 0.05
#' )
#'
#' Z <- c("t", "subtype")
#' Z_table <- unique(example[c(Z, probabilities)])
#' robin_wt(
#'   formula = y ~ xb + xc,
#'   data = example[setdiff(names(example), probabilities)],
#'   treatment = "treatment",
#'   probabilities = probabilities,
#'   Z = Z,
#'   Z_table = Z_table,
#'   treatments_for_compare = c("trt.1", "trt.2"),
#'   contrast = "difference",
#'   contrast_jac = NULL,
#'   family = gaussian(),
#'   stabilize = TRUE,
#'   alpha = 0.05
#' )
robin_wt <- function(formula, data, treatment, probabilities, Z = NULL, Z_table = NULL, treatments_for_compare,
                     contrast = "difference", contrast_jac=NULL, family=gaussian(), stabilize=TRUE, alpha=0.05,...) {
  if(is.null(probabilities)) stop("Assignment probabilities MUST be provided.")
  assert_subset(treatments_for_compare, probabilities)

  robin_estimate(formula = formula, data = data, treatment = treatment,
                 probabilities = probabilities, Z = Z, Z_table = Z_table,
                 treatments_for_compare = treatments_for_compare,
                 contrast = contrast, contrast_jac = contrast_jac, stabilize = stabilize,
                 alpha = alpha, method = "wt", post_strata = NULL, ...)
}

#' Post-Stratification Based Inference
#'
#' Provides robust inference methods via post stratification.
#'
#' @inheritParams robin_wt
#' @param post_strata (`character`) A string name of post-stratification variable. Default: `NULL`
#' @details
#' If family is `MASS::negative.binomial(NA)`, the function will use `MASS::glm.nb` instead of `glm`.
#' @export
#' @examples
#' probabilities <- c("trt.1", "trt.2", "trt.3", "trt.4")
#' robin_ps(
#'   formula = y ~ xb + xc,
#'   data = example,
#'   treatment = "treatment",
#'   probabilities = probabilities,
#'   Z = NULL,
#'   Z_table = NULL,
#'   post_strata = NULL,
#'   treatments_for_compare = c("trt.1", "trt.2"),
#'   contrast = "difference",
#'   contrast_jac = NULL,
#'   family = gaussian(),
#'   alpha = 0.05
#' )
#'
#' Z <- c("t", "subtype")
#' Z_table <- unique(example[c(Z, probabilities)])
#' robin_ps(
#'   formula = y ~ xb + xc,
#'   data = example[setdiff(names(example), probabilities)],
#'   treatment = "treatment",
#'   probabilities = c("trt.1", "trt.2", "trt.3", "trt.4"),
#'   Z = Z,
#'   Z_table = Z_table,
#'   post_strata = NULL,
#'   treatments_for_compare = c("trt.1", "trt.2"),
#'   contrast = "difference",
#'   contrast_jac = NULL,
#'   family = gaussian(),
#'   alpha = 0.05
#' )
robin_ps <- function(formula, data, treatment, probabilities = NULL, Z = NULL, Z_table = NULL, post_strata = NULL, treatments_for_compare,
                     contrast = "difference", contrast_jac = NULL, family = gaussian(), alpha = 0.05, ...) {

  robin_estimate(formula = formula, data = data, treatment = treatment,
                 probabilities = probabilities, Z = Z, Z_table = Z_table,
                 post_strata = post_strata, treatments_for_compare = treatments_for_compare,
                 contrast = contrast, contrast_jac = contrast_jac, family=family,
                 alpha=alpha, method = "ps", ...)
}
