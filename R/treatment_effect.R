#' Treatment Effect
#' @description Obtain treatment effect and variance from counter-factual prediction
#'
#' @param object Object from which to obtain treatment effect. Must be obtained from `estimate_effect()`.
#' @param pair (`integer` or `character`) Names or index of the treatment levels.
#' @param eff_measure (`function`) Treatment effect measurement function.
#' @param eff_jacobian (`function`) Treatment effect jacobian function.
#' @param alpha Nominal level
#' @param ... Additional arguments passed to `glm`
#'
#' @return A list of `treatment_effect` object with following elements:
#' - `mm_name`: name of the treatments to compare.
#' - `marginal_mean`: estimate of the treatment effect.
#' - `mmvariance`: estimate of the covariance matrix.
#' - `trt_effect`: estimate of the contrast.
#' - `variance`: estimate of the variance of contrast.
#' - `contrast`: name of the contrast function.
#' - `settings`: estimation settings.
#'
#' @export
treatment_effect <- function(object, pair, eff_measure, eff_jacobian, alpha,...) {
  UseMethod("treatment_effect", object)
}

#' @export
treatment_effect.prediction_cf <- function(
    object, pair = names(object), eff_measure, eff_jacobian, alpha,...) {

  checkmate::assert_function(eff_measure)
  contrast_name = deparse(substitute(eff_measure))
  contrast = ifelse(contrast_name=="h_diff", "difference",
                    ifelse(contrast_name=="h_ratio", "ratio",
                           ifelse(contrast_name=="h_odds_ratio", "odds ratio", "customized")))
  estimate <- object$estimation
  if (missing(pair)) {
    pair <- names(estimate$estimate)
  }
  checkmate::assert_vector(pair)
  checkmate::assert(
    checkmate::test_subset(pair, names(estimate$estimate)),
    checkmate::test_integerish(pair, lower = 1L, upper = length(estimate$estimate))
  )
  if (checkmate::test_integerish(pair)) {
    pair <- names(estimate$estimate)[pair]
  }
  trt_effect <- unname(eff_measure(estimate$estimate[pair]))

  inner_variance <- estimate$inner_variance
  if (missing(eff_jacobian)) {
    trt_jac <- numDeriv::jacobian(eff_measure, estimate[pair])
  } else {
    assert_function(eff_jacobian)
    trt_jac <- eff_jacobian(estimate$estimate[pair])
  }
  trt_var <- trt_jac %*% inner_variance %*% t(trt_jac)


  pair_names <- outer(pair, pair, FUN = paste, sep = " - ")
  structure(
    list(
      mm_name = pair,
      marginal_mean = estimate$estimate,
      mmvariance = estimate$inner_variance,
      trt_effect = trt_effect,
      variance = diag(trt_var),
      contrast = contrast,
      settings = attr(object, "settings")
    ),
    name = pair_names[lower.tri(pair_names)],
    fit.j = object$fit.j,
    treatment = attr(object, "treatment_name"),
    prob_mat = attr(object, "prob_mat"),
    Z = attr(object, "Z"),
    sample_size = object$sample_size,
    alpha = alpha,
    post_strata = attr(object, "post_strata"),
    data = attr(object, "data"),
    method = estimate$method,
    class = "treatment_effect"
  )
}

#' @rdname treatment_effect
difference <- function(object, ...) {
  treatment_effect(object, eff_measure = h_diff, eff_jacobian = h_jac_diff, ...)
}
#' @rdname treatment_effect
risk_ratio <- function(object, ...) {
  treatment_effect(object, eff_measure = h_ratio, eff_jacobian = h_jac_ratio, ...)
}
#' @rdname treatment_effect
odds_ratio <- function(object, ...) {
  treatment_effect(object, eff_measure = h_odds_ratio, eff_jacobian = h_jac_odds_ratio, ...)
}

#' Contrast Functions and Jacobians
#' @rdname contrast
#' @param x (`numeric`) Vector of values.
#' @return Vector of contrasts, or matrix of jacobians.
#' @examples
#' h_diff(1:3)
#' h_jac_ratio(1:3)
#' @export
h_diff <- function(x) {
  assert_numeric(x)
  d <- outer(x, x, `-`)
  d[lower.tri(d)]
}

#' @rdname contrast
#' @export
h_jac_diff <- function(x) {
  assert_numeric(x)
  n <- length(x)
  l <- h_lower_tri_idx(n)
  ret <- matrix(0, nrow = nrow(l), ncol = n)
  ret[cbind(seq_len(nrow(ret)), l[, 1])] <- 1
  ret[cbind(seq_len(nrow(ret)), l[, 2])] <- -1
  ret
}

#' @rdname contrast
#' @export
h_ratio <- function(x) {
  assert_numeric(x, lower = 0)
  d <- outer(x, x, `/`)
  d[lower.tri(d)]
}

#' @rdname contrast
#' @export
h_jac_ratio <- function(x) {
  assert_numeric(x, lower = 0)
  n <- length(x)
  l <- h_lower_tri_idx(n)
  ret <- matrix(0, nrow = nrow(l), ncol = n)
  ret[cbind(seq_len(nrow(ret)), l[, 1])] <- 1 / x[l[, 2]]
  ret[cbind(seq_len(nrow(ret)), l[, 2])] <- -x[l[, 1]] / x[l[, 2]]^2
  ret
}

#' @rdname contrast
#' @export
h_odds_ratio <- function(x) {
  assert_numeric(x, lower = 0, upper = 1)
  y <- x / (1 - x)
  h_ratio(y)
}

#' @rdname contrast
#' @export
h_jac_odds_ratio <- function(x) {
  assert_numeric(x, lower = 0)
  n <- length(x)
  l <- h_lower_tri_idx(n)
  ret <- matrix(0, nrow = nrow(l), ncol = n)
  ret[cbind(seq_len(nrow(ret)), l[, 1])] <- (1 - x[l[, 2]]) / ((1 - x[l[, 1]])^2 * x[l[, 2]])
  ret[cbind(seq_len(nrow(ret)), l[, 2])] <- -x[l[, 1]] / ((1 - x[l[, 1]]) * x[l[, 2]]^2)
  ret
}

#' Lower Triangular Index
#' @param n (`int`) Number of rows/columns.
#' @return Matrix of lower triangular indices.
#' @keywords internal
h_lower_tri_idx <- function(n) {
  rc <- c(n, n)
  which(.row(rc) > .col(rc), arr.ind = TRUE)
}

#' @export
print.treatment_effect <- function(x, ...) {
  alpha <- attr(x,"alpha")
  data <- attr(x, "data")
  post_strata <- attr(x, "post_strata")
  prob_mat <- round(attr(x, "prob_mat"), 2)
  settings <- x$settings
  stratify_by <- settings$stratify_by
  Z <- attr(x, "Z")


  cat("Method: ", attr(x, "method"), "\n")
  if(settings$method=="ps"){
    if(!is.null(stratify_by)) {
      cat(paste0("Post stratification is done by variable ", stratify_by, " specified by stratify_by.\n"))
          } else {cat(paste0("Post stratification is done by the joint levels of the randomization variables specified by randomization_var_colnames.\n"))}
  }
  cat("Model : ", deparse(as.formula(attr(x, "fit"))), "\n")
  cat("Family: ", attr(x, "fit")$family[[1]], "\n")

  if(settings$estimated_propensity) {cat("Estimated Propensity Score is used.\n")}
  p <- ""
  if(settings$method == "wt" & settings$estimated_propensity) {
    p <- "Estimated "
  } else if(settings$method=="ps" & (!is.null(stratify_by) | is.null(stratify_by) & !settings$rand_table)) p="Estimated "
  cat(paste0(p,"Randomization Probabilities (among the entire concurrent and eligible (ECE) population):"), "\n")
  prob_tab <- prob_table_generate(data, post_strata, prob_mat, Z)
  row.names(prob_tab) <- 1:nrow(prob_tab)

  prob_tab_output_lines <- capture.output(print.data.frame(prob_tab, row.names = TRUE)) # Force row.names printing
  cat(paste(prob_tab_output_lines, collapse = "\n"), "\n")

  cat("\n")
  cat("Nominal Level: ", alpha, "\n")
  cat("---------------------------\n")
  cat("Marginal Mean: \n")

  mm <- as.numeric(x$marginal_mean)
  trt_sd <- sqrt(diag(x$mmvariance))
  z_value <- mm / trt_sd
  ci.lb <- mm - qnorm(1-alpha/2) * trt_sd
  ci.ub <- mm + qnorm(1-alpha/2) * trt_sd
  coef_mat <- matrix(
    c(
      mm,
      trt_sd,
      z_value,
      ci.lb,
      ci.ub
    ),
    nrow = length(mm)
  )
  colnames(coef_mat) <- c("Estimate", "Std.Err","Z Value", "lower.CL", "upper.CL")
  row.names(coef_mat) <- x$mm_name
  stats::printCoefmat(
    coef_mat,
    cs.ind = c(1,2,4,5),
    tst.ind = 3,
    digits = 3,
    P.values = FALSE
  )
  cat("\n")
  cat("---------------------------\n")
  cat("Treatment Effect: \n")
  cat("Contrast: ", x$contrast,"\n")

  trt_sd <- sqrt(x$variance)
  trt_effect <- x$trt_effect
  z_value <- as.numeric(trt_effect) / trt_sd
  ci.lb <- trt_effect - qnorm(1-alpha/2) * trt_sd
  ci.ub <- trt_effect + qnorm(1-alpha/2) * trt_sd
  p <- 2 * pnorm(abs(z_value), lower.tail = FALSE)
  coef_mat <- matrix(
    c(
      trt_effect,
      trt_sd,
      z_value,
      ci.lb,
      ci.ub,
      p
    ),
    nrow = length(trt_effect)
  )
  colnames(coef_mat) <- c("Estimate", "Std.Err", "Z Value","lower.CL", "upper.CL",  "Pr(>|z|)")
  row.names(coef_mat) <- attr(x, "name")
  stats::printCoefmat(
    coef_mat,
    cs.ind = c(1,2,4,5),
    tst.ind = 3,
    digits = 3,
    P.values = TRUE
  )
}
