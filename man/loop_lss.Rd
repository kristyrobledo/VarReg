% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loop_lss.R
\name{loop_lss}
\alias{loop_lss}
\title{The EM loop for the LSS model}
\usage{
loop_lss(
  alldat,
  xiold,
  omega2old,
  nuold,
  mean.ind,
  var.ind,
  nu.ind,
  para.space,
  maxit,
  eps,
  int.maxit,
  print.it
)
}
\arguments{
\item{alldat}{Dataframe containing all the data for the models. Outcome in the first column.}

\item{xiold}{Vector of initial location parameter estimates to be fit in the location model.}

\item{omega2old}{Vector of initial scale2 parameter estimates to be fit in the scale2 model.}

\item{nuold}{Vector of initial nu parameter estimates to be fit in the nu model.}

\item{mean.ind}{Vector containing the column numbers of the data in 'alldat' to be fit as covariates in the location model.}

\item{var.ind}{Vector containing the column numbers of the data in 'alldat' to be fit as covariates in the scale2 model. FALSE indicates a constant variance model.}

\item{nu.ind}{Vector containing the column numbers of the data in 'alldat' to be fit as covariates in the nu model. NULL indicates constant model.}

\item{para.space}{Parameter space to search for variance parameter estimates. "positive" means only search positive parameter space, "negative" means search only negative parameter space and "all" means search all.}

\item{maxit}{Number of maximum iterations for the main EM algorithm.}

\item{eps}{Very small number for the convergence criteria.}

\item{int.maxit}{Number of maximum iterations for the internal EM algorithm for the location and scale.}

\item{print.it}{Logical to indicate if the estimates for each iteration should be printed.}
}
\value{
A list of the results from the algorithm, including conv, reldiff, it, mean, xi.new, omega2.new, nu.new, fitted.xi
\itemize{
\item\code{conv}: Logical argument indicating if convergence occurred
\item\code{it}: Total iterations performed of the EM algorithm
 \item\code{reldiff}: the positive convergence tolerance that occured at the final iteration
 \item\code{xinew}: Vector of location parameter estimates
 \item\code{omega2new}: Vector of scale squared parameter estimates
 \item\code{nunew}: Vector of shape parameter estimates
 \item\code{fitted.xi}: Vector of fitted location estimates
 }
}
\description{
\code{loop_lss} is the EM loop function for the LSS model to be utilised by various other higher level functions
}
