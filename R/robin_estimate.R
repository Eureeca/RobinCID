#' @noRd
robin_estimate <- function(data,
                           estimand = list(tx_colname = NULL,
                                           tx_to_compare = NULL),
                           design = list(randomization_var_colnames = NULL,
                                         randomization_table = NULL),
                           stratify_by, estimated_propensity,
                           outcome_model = list(formula = NULL,
                                                family = gaussian()),
                           contrast_specs = list(contrast = "difference",
                                                 contrast_jac = NULL),
                           alpha, method,
                           ...){
  validate_inputs(estimand, design, outcome_model, estimated_propensity, method)
  treatment <- estimand$tx_colname
  treatments_for_compare <- estimand$tx_to_compare

  post_strata <- stratify_by

  formula <-outcome_model$formula
  family <- outcome_model$family

  contrast <- contrast_specs$contrast
  contrast_jac <- contrast_specs$contrast_jac

  treatment_names <- unique(data[[treatment]])

  assert_subset(all.vars(formula), names(data))
  assert_subset(treatment, names(data))
  assert_subset(as.character(treatments_for_compare), as.character(data[[treatment]]))

  Z <- design$randomization_var_colnames

  consistency_check(data, estimand, design, stratify_by)

  data <- assign_prob_and_strata(
    data = data,
    estimand = estimand,
    design = design,
    method = method,
    estimated_propensity = estimated_propensity,
    stratify_by = stratify_by
  )
  prob_mat <- data[treatments_for_compare]


  data[[treatment]] <- as.factor(data[[treatment]])
  treatments_for_compare <- factor(treatments_for_compare,
                                   levels = levels(data[[treatment]]),
                                   labels = levels(data[[treatment]]))

  # ECE sample
  ECE_subset <- apply(data[as.character(treatments_for_compare)] > 0, 1, all)
  data <- data[ECE_subset, ]
  if (nrow(data) < 1)
    stop("ECE sample size is 0!")
  prob_mat <- prob_mat[ECE_subset, ]
  attr(prob_mat, "Z") <- Z

  if (identical(family$family, "Negative Binomial(NA)")) {
    fit.j <-  MASS::glm.nb(formula, family = family, data = data[data[[treatment]]==treatments_for_compare[1],], ...)
    fit.k <-  MASS::glm.nb(formula, family = family, data = data[data[[treatment]]==treatments_for_compare[2],], ...)
  } else {
    fit.j <-  glm(formula, family = family, data = data[data[[treatment]]==treatments_for_compare[1],], ...)
    fit.k <-  glm(formula, family = family, data = data[data[[treatment]]==treatments_for_compare[2],], ...)
  }
  settings <- list(method = method,
                  estimated_propensity = estimated_propensity,
                  stratify_by = stratify_by,
                  rand_table = !is.null(design$randomization_table))
  pc <- predict_counterfactual(fit.j = fit.j, fit.k = fit.k, treatment = treatment,
                               treatments_for_compare = treatments_for_compare,
                               prob_mat = prob_mat, post_strata = post_strata,
                               data = data, stabilize = TRUE, settings = settings)

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
#' Provides robust inference via inverse probability weighting.
#'
#' @param data (`data.frame`) A data frame containing the dataset.
#' @param estimand (`list`) A list specifying the estimand, with two elements:
#'   - `tx_colname` (`character`): The column name of the treatment variable in `data`.
#'   - `tx_to_compare` (`character vector`): A vector specifying exactly two treatment levels to compare.
#' @param design (`list`) A list specifying randomization information, with two elements:
#'   - `randomization_var_colnames` (`character vector`): Column names of randomization variables in `data`.
#'   - `randomization_table` (`data.frame`, default: `NULL`): A data frame containing treatment assignment probabilities
#'     for each level of the randomization variables. See *Details*.
#' @param estimated_propensity (`logical`, default: `TRUE`) Whether to use estimated propensity scores.
#' @param outcome_model (`list`) A list specifying the outcome working model, with two elements:
#'   - `formula` (`formula`): The regression formula for the analysis.
#'   - `family` A description of the error distribution and link function for the model. Default: `gaussian()`.
#' @param contrast_specs (`list`) A list specifying the contrast function and its Jacobian:
#'   - `contrast` (`function` or `character`): A function to compute the treatment effect, or one of `"difference"`,
#'     `"risk_ratio"`, or `"odds_ratio"` for default contrasts.
#'   - `contrast_jac` (`function`, optional): A function to compute the Jacobian of the contrast function.
#'     Ignored if using default contrasts.
#' @param alpha (`numeric`) The nominal significance level. Default: `0.05`.
#' @param ... Additional arguments passed to `glm`.
#' @details
#' If `randomization_table` is provided, it must include columns corresponding to `randomization_var_colnames`,
#' as well as treatment assignment probability columns named after the treatment levels in `tx_colname` from `data`.
#'
#' If `family` is `MASS::negative.binomial(NA)`, the function will use `MASS::glm.nb` instead of `glm`.
#'
#' @export
#'
#' @return A treatment_effect object.
#' @examples
#' data_sim <- RobinCID::example
#' tx_colname <- "treatment"
#' treatment_levels <- unique(data_sim[[tx_colname]])
#' tx_to_compare <- c("trt.1", "trt.3")
#' randomization_var_colnames <- c("t", "subtype")
#' df <- data_sim[c("xb", "xc", tx_colname, randomization_var_colnames, "y")]
#' randomization_table <- unique(data_sim[c(randomization_var_colnames, treatment_levels)])
#' robin_wt(
#'   data = df,
#'   estimand = list(tx_colname = tx_colname,
#'                   tx_to_compare = tx_to_compare),
#'   design = list(randomization_var_colnames = randomization_var_colnames,
#'                 randomization_table = randomization_table),
#'   estimated_propensity = FALSE,
#'   outcome_model = list(formula = y ~ 1,
#'                        family = gaussian())
#' )
#'
robin_wt <- function(data,
                     estimand = list(tx_colname = NULL,
                                     tx_to_compare = NULL),
                     design = list(randomization_var_colnames = NULL,
                                   randomization_table = NULL),
                     estimated_propensity = TRUE,
                     outcome_model = list(formula = NULL,
                                          family = gaussian()),
                     contrast_specs = list(contrast = "difference",
                                           contrast_jac = NULL),
                     alpha=0.05,
                     ...) {

  robin_estimate(data = data,
                 estimand = estimand,
                 design = design,
                 estimated_propensity = estimated_propensity,
                 outcome_model = outcome_model,
                 contrast_specs = contrast_specs,
                 alpha=alpha, method = "wt", stratify_by = NULL, ...)
}

#' Post-Stratification Based Inference
#'
#' Provides robust inference via post stratification.
#'
#' @inheritParams robin_wt
#' @param stratify_by (`character`, optional) The column name of the stratification variable in `data`.
#' If provided, `stratify_by` overrides `design`.
#' @details
#' If family is `MASS::negative.binomial(NA)`, the function will use `MASS::glm.nb` instead of `glm`.
#' @export
#' @return A treatment_effect object.
#' @examples
#' data_sim <- RobinCID::example
#' tx_colname <- "treatment"
#' treatment_levels <- unique(data_sim[[tx_colname]])
#' tx_to_compare <- c("trt.1", "trt.3")
#' randomization_var_colnames <- c("t", "subtype")
#' df <- data_sim[c("xb", "xc", tx_colname, randomization_var_colnames, "y")]
#' randomization_table <- unique(data_sim[c(randomization_var_colnames, treatment_levels)])
#' robin_ps(
#'   data = df,
#'   estimand = list(tx_colname = tx_colname,
#'                   tx_to_compare = tx_to_compare),
#'   design = list(randomization_var_colnames = randomization_var_colnames,
#'                 randomization_table = randomization_table),
#'   stratify_by = NULL,
#'   outcome_model = list(formula = y ~ 1,
#'                        family = gaussian())
#' )
robin_ps <- function(data,
                     estimand = list(tx_colname = NULL,
                                     tx_to_compare = NULL),
                     design = list(randomization_var_colnames = NULL,
                                   randomization_table = NULL),
                     stratify_by = NULL,
                     outcome_model = list(formula = NULL,
                                          family = gaussian()),
                     contrast_specs = list(contrast = "difference",
                                           contrast_jac = NULL),
                     alpha=0.05,
                     ...) {

  robin_estimate(data = data,
                 estimand = estimand,
                 design = design,
                 stratify_by = stratify_by,
                 outcome_model = outcome_model,
                 contrast_specs = contrast_specs,
                 alpha=alpha, method = "ps", estimated_propensity = FALSE, ...)
}
