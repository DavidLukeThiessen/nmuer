% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nmuer.R
\name{nmuer}
\alias{nmuer}
\title{Non-Monotone Unified Estimate}
\usage{
nmuer(
  main_model,
  working_models,
  data,
  main_weight_model = NULL,
  working_weight_models = NULL,
  variance_estimate = c("bootstrap", "sandwich"),
  num_boots = 500,
  DEBUG = FALSE
)
}
\arguments{
\item{main_model}{(Required) This is the main analysis model. This should
be either a glm() or coxph() function which includes a formula but does
not include a data argument}

\item{working_models}{(Required) A list containing the working models,
similar to the main_model specification.}

\item{data}{(Required) A dataframe with 1 row for each individual.}

\item{main_weight_model}{(Optional) A function that takes at minimum a
Data argument and a observation_vector argument, and outputs a vector of
weights. These weights are used in the main_model.}

\item{working_weight_models}{(Optional) A list of functions that each
take DATA and observation_vector arguments and output a vector of
weights.}

\item{variance_estimate}{(Optional) There are two different ways of
estimating the variance. Only "bootstrap" is currently allowed. 
"bootstrap" uses bootstrap resampling on rows of the data to estimate
the variance. "sandwich" will use a robust sandwich estimate of the form
A^-1 * B * A^T-1, where A is the matrix of partial
derivatives of the score with respect to the parameters and B is the
matrix of crossproducts of the score vectors.}

\item{num_boots}{(Optional) The number of bootstrap replications to use
when using the bootstrap method of estimating variance. 200 was
recommended as the minimum by Efron, pg 52.}

\item{DEBUG}{(Optional) Prints extra model fitting and call information}
}
\description{
This function fits the unified estimate for either a Cox
or a GLM model. It uses bootstrapping to estimate variance.
}
\details{
TODO
}
\examples{
if(require(mice)) {
  nmuer_glm_mod <- nmuer(
    main_model = glm(hm ~ age*sex),
    working_models = list(glm(hr ~ age*sex)),
    data = mice::selfreport)
  summary(nmuer_glm_mod)
}
}
