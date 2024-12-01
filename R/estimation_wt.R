#' Compute Estimates and Covariance Matrix
#'
#' @param ret counterfactual prediction
#' @param treatments_for_compare description
#' @param data (`data.frame`) data
#' @param prob_mat (`data.frame`) treatment assignment probability
#' @param stabilize (`logical`) whether to stabilize
#' @param y Observed outcome
#' @param treatment name of treatment
#'
#' @export
estimation_wt <- function(ret, y, treatment, treatments_for_compare, data, prob_mat, stabilize=T){

  pij <- prob_mat[[treatments_for_compare[1]]]
  pik <- prob_mat[[treatments_for_compare[2]]]

  njk <- nrow(data)

  A.j <- as.numeric(data[[treatment]] == treatments_for_compare[1])
  A.k <- as.numeric(data[[treatment]] == treatments_for_compare[2])

  weightj <- if (stabilize) sum(A.j / pij) else njk
  weightk <- if (stabilize) sum(A.k / pik) else njk

  pred.jk = ret[, 1]
  pred.kj = ret[, 2]

  theta.jk <- sum(A.j * (y - pred.jk) / pij) / weightj + sum(pred.jk) / njk
  theta.kj <- sum(A.k * (y - pred.kj) / pik) / weightk + sum(pred.kj) / njk

  delta.jk = sum(A.j * (y-pred.jk) / pij) / njk
  delta.kj = sum(A.k * (y-pred.kj) / pik) / njk

  q.jk.j = cov((y-pred.jk)[A.j==1], pred.jk[A.j==1])
  q.kj.k = cov((y-pred.kj)[A.k==1], pred.kj[A.k==1])
  q.jk.k = cov((y-pred.jk)[A.k==1], pred.jk[A.k==1])
  q.kj.j = cov((y-pred.kj)[A.j==1], pred.kj[A.j==1])
  q.jk = cov(pred.jk, pred.kj)

  sigma.jk = var(pred.jk)
  sigma.kj = var(pred.kj)

  lambda.jk = 2 * q.jk.j + sigma.jk
  lambda.kj = 2 * q.kj.k + sigma.kj

  c.jk = q.kj.j + q.jk.k + q.jk

  sigma.11 = sum(A.j*(y-pred.jk - delta.jk)^2/pij^2)/njk + lambda.jk
  sigma.22 = sum(A.k*(y-pred.kj - delta.kj)^2/pik^2)/njk + lambda.kj
  sigma.12 = c.jk

  inner_variance = matrix(c(sigma.11, sigma.12,sigma.12,sigma.22),ncol = 2, nrow = 2) / njk
  estimate = c(theta.jk, theta.kj)
  names(estimate) = treatments_for_compare

  list(estimate = estimate, inner_variance = inner_variance)
}
