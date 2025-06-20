% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{assign_prob_and_strata}
\alias{assign_prob_and_strata}
\title{Assign Probability according to Design}
\usage{
assign_prob_and_strata(
  data,
  estimand,
  design = list(randomization_var_colnames = NULL, randomization_table = NULL),
  method,
  estimated_propensity = TRUE,
  stratify_by = NULL
)
}
\arguments{
\item{data}{(\code{data.frame}) Input data frame.}

\item{estimand}{(\code{list}) A list specifying the estimand.}

\item{design}{(\code{list}) A list describing the randomization design. See \code{Details}.}

\item{method}{estimation method.}

\item{estimated_propensity}{Whether to use estimated propensity score.}

\item{stratify_by}{The column name of stratification variable in \code{data}.}
}
\value{
A new \code{data} with columns of the treatment assignment probability.
}
\description{
Assign Probability according to Design
}
\details{
\code{design} has two elements: \code{randomization_var_colnames} (\code{vector}) and \code{randomization_table} (\code{data.frame})
}
