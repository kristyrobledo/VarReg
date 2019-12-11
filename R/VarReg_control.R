#' Auxillary for controlling VarReg fitting
#'
#' Use \code{VarReg.control} to determine parameters for the fitting of \code{\link{semiVarReg}}. Typically only used internally within functions.
#' @param bound.tol Positive tolerance for specifying the interior of the parameter space. This allows the algorithm to terminate early if an interior maximum is found. If set to \code{bound.tol=Inf}, no early termination is attempted.
#' @param epsilon  Positive convergence tolerance. If \eqn{\theta} is a vector of estimates, convergence is declared when \eqn{\sqrt{(\sum (\theta_{old} - \theta_{new})^2)}/ \sqrt{\sum (\theta_{old})^2} }. This should be smaller than \code{bound.tol}.
#' @param maxit integer giving the maximum number of EM algorithm iterations for a given parameterisation.
#' @return A list of the three components: \code{bound.tol}, \code{epsilon} and \code{maxit} .
#'
#'@details This is used similarly to \code{\link[stats]{glm.control}}. If required, it may be internally passed to another function.

#' @export

VarReg.control<-function (bound.tol = 1e-05, epsilon = 1e-06, maxit = 1000) {
  if (!is.numeric(bound.tol) || bound.tol <= 0)
    stop("value of 'bound.tol' must be > 0")
  if (!is.numeric(epsilon) || epsilon <= 0)
    stop("value of 'epsilon' must be > 0")
  if (epsilon >= bound.tol)
    warning("'epsilon' should be smaller than 'bound.tol'",
            call. = FALSE)
  if (!is.numeric(maxit) || maxit <= 0)
    stop("maximum number of iterations must be > 0")
  list(bound.tol = bound.tol, epsilon = epsilon, maxit = maxit)
}