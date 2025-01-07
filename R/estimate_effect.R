#' Compute Estimates and Covariance Matrix
#'
#' @param ret counterfactual prediction
#' @param treatments_for_compare description
#' @param data (`data.frame`) data
#' @param prob_mat (`data.frame`) treatment assignment probability
#' @param post_strat (`character`) A string name of post-stratification variable
#' @param stabilize (`logical`) whether to stabilize
#' @param y Observed outcome
#' @param treatment name of treatment
#'
estimate_effect <- function(ret, y, treatment, treatments_for_compare, data, prob_mat, post_strat, stabilize) {
  UseMethod("estimate_effect", ret)
}

#' @export
estimate_effect.wt <- function(ret, y, treatment, treatments_for_compare, data, prob_mat, post_strat, stabilize){

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

  list(estimate = estimate, inner_variance = inner_variance,
       method="Inverse Probability Weighting")
}

#' @export
estimate_effect.ps <- function(ret, y, treatment, treatments_for_compare, data, prob_mat, post_strat, stabilize) {


  njk <- nrow(data)

  A.j <- (data[[treatment]] == treatments_for_compare[1])
  A.k <- (data[[treatment]] == treatments_for_compare[2])

  pred.jk <- ret[, 1]
  pred.kj <- ret[, 2]

  pij <- prob_mat[[treatments_for_compare[1]]]
  pik <- prob_mat[[treatments_for_compare[2]]]

  temp_df <- data.frame(pij = pij, pik = pik, A.j = A.j, A.k = A.k, pred.j = pred.jk, pred.k = pred.kj, y = y)

  if(is.null(post_strat)){
    unique_pairs <- unique(temp_df[, c("pij", "pik")])
    n_pairs <- nrow(unique_pairs)
    strata <- 1:n_pairs
    temp_df$post_strat <- match(interaction(temp_df$pij, temp_df$pik), interaction(unique_pairs$pij, unique_pairs$pik))
  } else {
    n_pairs <- length(unique(data[[post_strat]]))
    temp_df$post_strat <- data[[post_strat]]
    strata <- unique(data[[post_strat]])
  }
  thetas.j <- thetas.k <- sigma11 <- sigma22 <- sigma12 <- Y.bar.j <- Y.bar.k <- ns <- numeric(n_pairs)

  for (i in seq_len(n_pairs)) {
    pair_data <- temp_df[temp_df$post_strat==strata[i], ]

    n.jk <- nrow(pair_data)

    theta.j <- sum(pair_data$A.j * (pair_data$y - pair_data$pred.j)) / sum(pair_data$A.j) + mean(pair_data$pred.j)
    theta.k <- sum(pair_data$A.k * (pair_data$y - pair_data$pred.k)) / sum(pair_data$A.k) + mean(pair_data$pred.k)

    y.bar.j <- mean(pair_data$y[pair_data$A.j])
    y.bar.k <- mean(pair_data$y[pair_data$A.k])

    pij.hat <- mean(pair_data$A.j)
    pik.hat <- mean(pair_data$A.k)

    sigma.jk <- var(pair_data$pred.j)
    sigma.kj <- var(pair_data$pred.k)

    q.jk.j <- cov((pair_data$y - pair_data$pred.j)[pair_data$A.j], pair_data$pred.j[pair_data$A.j])
    q.kj.k <- cov((pair_data$y - pair_data$pred.k)[pair_data$A.k], pair_data$pred.k[pair_data$A.k])
    q.jk.k <- cov((pair_data$y - pair_data$pred.j)[pair_data$A.k], pair_data$pred.j[pair_data$A.k])
    q.kj.j <- cov((pair_data$y - pair_data$pred.k)[pair_data$A.j], pair_data$pred.k[pair_data$A.j])
    q.jk <- cov(pair_data$pred.j, pair_data$pred.k)

    lambda.jk <- 2 * q.jk.j + sigma.jk
    lambda.kj <- 2 * q.kj.k + sigma.kj
    c.jk <- q.kj.j + q.jk.k + q.jk
    tau.jk <- var(pair_data$y[pair_data$A.j] - pair_data$pred.j[pair_data$A.j])
    tau.kj <- var(pair_data$y[pair_data$A.k] - pair_data$pred.k[pair_data$A.k])
    sigma.11 <- tau.jk / pij.hat + lambda.jk
    sigma.22 <- tau.kj / pik.hat + lambda.kj
    sigma.12 <- c.jk
    thetas.j[i] <- theta.j
    thetas.k[i] <- theta.k
    sigma11[i] <- sigma.11
    sigma22[i] <- sigma.22
    sigma12[i] <- sigma.12
    Y.bar.j[i] <- y.bar.j
    Y.bar.k[i] <- y.bar.k
    ns[i] <- n.jk
  }

  theta.jk <- weighted.mean(x = thetas.j, w = ns)
  theta.kj <- weighted.mean(x = thetas.k, w = ns)
  Y.bar.js <- rep(Y.bar.j, ns)
  Y.bar.ks <- rep(Y.bar.k, ns)
  Gamma.jk <- cov(cbind(Y.bar.js, Y.bar.ks))
  Sigma.11 <- weighted.mean(x = sigma11, w = ns)
  Sigma.12 <- weighted.mean(x = sigma12, w = ns)
  Sigma.22 <- weighted.mean(x = sigma22, w = ns)
  Sigma.jk <- matrix(c(
    Sigma.11, Sigma.12,
    Sigma.12, Sigma.22
  ), 2, 2) + Gamma.jk

  inner_variance <- Sigma.jk / sum(ns)
  estimate <- c(theta.jk, theta.kj)
  names(estimate) <- treatments_for_compare
  list(estimate = estimate, inner_variance = inner_variance,
       method = "Post Stratification")
}


