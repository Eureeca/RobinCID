#' Title
#'
#' @param ret counterfactual prediction
#' @param treatments_for_compare description
#' @param data (`data.frame`) data
#' @param prob_mat (`data.frame`) treatment assignment probability
#' @param stabilize (`logical`) whether to stabilize
#' @param y Observed outcome
#' @param treatment name of treatment
#' @param ...
#'
#' @export
estimation_wt <- function(ret, y, treatment, treatments_for_compare, data, prob_mat, stabilize=T,...){

  pij <- prob_mat[[treatments_for_compare[1]]]
  pik <- prob_mat[[treatments_for_compare[2]]]

  njk <- nrow(data)

  A.j = data[[treatment]] == treatments_for_compare[1]
  A.k = data[[treatment]] == treatments_for_compare[2]

  weightj = ifelse(stabilize, sum(A.j / pij), njk)
  weightk = ifelse(stabilize, sum(A.k / pik), njk)

  pred.jk = ret[,treatments_for_compare[1]]
  pred.kj = ret[,treatments_for_compare[2]]

  theta.jk <- sum(A.j * (y - pred.jk) / pij) / weightj + sum(pred.jk) / njk
  theta.kj <- sum(A.k * (y - pred.kj) / pik) / weightk + sum(pred.kj) / njk

  c(theta.jk, theta.kj)
}
