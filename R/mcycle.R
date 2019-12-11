#' mcycle dataset.
#'
#' A dataset containing 133 observations from a simulated motorcycle accident, used to test crash helmets.
#'
#' @format A data frame with 133 rows and 2 variables:
#' \describe{
#'   \item{times}{in milliseconds from time of impact}
#'   \item{accel}{in g, acceleration of the head}
#'   ...
#' }
#' @source Silverman, B. W. (1985) Some aspects of the spline smoothing approach to non-parametric curve
#'fitting. Journal of the Royal Statistical Society series B 47, 1-52.


#' @references Venables, W. N. and Ripley, B. D. (1999) Modern Applied Statistics with S-PLUS. Third Edition.
#' Springer.
#'@examples
#'library(VarReg)
#'data(mcycle)
#'attach(mcycle)
#'plot(times,accel)
"mcycle"