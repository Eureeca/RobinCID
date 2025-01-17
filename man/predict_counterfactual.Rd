% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_counterfactual.R
\name{predict_counterfactual}
\alias{predict_counterfactual}
\title{Counterfactual Prediction}
\usage{
predict_counterfactual(
  fit.j,
  fit.k,
  treatment,
  treatments_for_compare,
  prob_mat,
  post_strat,
  data,
  stabilize,
  method
)
}
\arguments{
\item{fit.j}{fitted object for trt j.}

\item{fit.k}{fitted object for trt k.}

\item{treatment}{name of treatment column}

\item{treatments_for_compare}{(\code{character}) Treatments for comparison}

\item{prob_mat}{(\code{data.frame}) treatment assignment probabilities}

\item{post_strat}{(\code{character}) A string name of post-stratification variable.}

\item{data}{(\code{data.frame}) raw dataset.}

\item{stabilize}{stabilize}

\item{method}{estimation method}
}
\value{
Numeric matrix of counter factual prediction.
}
\description{
Obtain counterfactual prediction of a fit.
}
