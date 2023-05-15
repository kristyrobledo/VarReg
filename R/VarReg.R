#' VarReg: Semi-parametric mean and variance regression
#'
#' Methods for fitting semi-parametric mean and variance models, with normal or censored data.
#' Also extended to allow a regression in the location, scale and shape parameters.
#'
#' This package provides functions to fit semi-parametric mean and variance regression models. These models
#' are based upon EM-type algorithms, which can have more stable convergence properties than other
#' algorithms for additive variance regression models.
#'
#'The primary function to use for linear and semi-parametric mean and variance models is \code{\link{semiVarReg}}.
#'This function also is able to fit models to censored outcome data. There is also a plot function for these
#'models called \code{\link{plotVarReg}}.
#'A search function has also been produced in order to assist users to find the optimal number of knots in
#'the model (\code{\link{searchVarReg}}).
#'
#' The other functions that are of particular use are \code{\link{lssVarReg}} and its plot function
#' \code{\link{plotlssVarReg}}. This uses the skew-normal distribution and combines the EM algorithm with
#' a coordinate-ascent type algorithm in order to fit a regression model in the location, scale and shape,
#'  therefore extending the semi-parametric models to non-normal data.
#'
#'   Multivariate models can be fit with \code{\link{semiVarReg.multi}} and \code{\link{lssVarReg.multi}}
#'
#'
#' @docType package
#' @name VarReg
#' @author Kristy Robledo \email{robledo.kristy@@gmail.com}
"_PACKAGE"