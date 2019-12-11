#' Plots graphics for a location, scale and shape regression model
#'
#' \code{plotlssVarReg} is used to produce graphics for models fit in the \code{VarReg} package with the
#'  \code{lssVarReg} function. As the skew-normal distribution is used to fit this type of model, the data needs
#'  to be transformed from the SN parameters (location, scale and shape) to the typical mean,
#'  variance and skew parameters.
#' @param x Object of class lssVarReg (output from \code{\link{lssVarReg}}).
#' @param knot.lines  Logical to show the knot lines on the graphics (if model is type "semi").
#' Default is \code{TRUE}
#' @param xlab Label to be placed on the x axis of graphics (covariate)
#' @param ylab Label to be placed on the y axis of graphics (outcome)
#' @return A graphic is returned, as well as a dataframe. The graphic returned is a 2 by 2 plot of:
#' \itemize{
#' \item the mean function over the x-variable, with or without the knot lines indicated
#' \item the variance function over the x-variable, with or without the knot lines indicated
#' \item the skew function over the x-variable, with or without the knot lines indicated
#' \item a Q-Q plot of the squared residuals from the model, plotted against the Chi-squared (df=1)
#'  distribution. For data from a skew-normal distribution, these residuals should follow a
#'  Chi-squared (df=1) distribution, regardless of skew.
#'  }
#'
#' The dataframe returned contains the following columns:
#' \itemize{
#' \item \code{x}: x variable
#' \item \code{y}: y variable
#' \item \code{eta}:  (\eqn{\eta}), the location parameter
#' \item \code{omega}: (\eqn{\omega}), the scale parameter
#' \item \code{shape}: (\eqn{\nu}),  the shape parameter
#' \item \code{predicted~mean}: (\eqn{\mu}), the mean
#' \item \code{predicted~variance}: (\eqn{\sigma^2}),  the variance
#' \item \code{predicted~skewness}: (\eqn{\gamma}),  the skew
#' \item \code{stand.res2}: the standardised residuals squared.
#' }
#'
#'@seealso \code{\link{lssVarReg}}
#'
#'@examples
#'data(mcycle)
#' ## not run. LSS model followed by the basic plot command
#' ##lssmodel<-lssVarReg(mcycle$accel, mcycle$times,  locationmodel="linear", scale2model="linear",
#' ##shapemodel="constant", maxit=10000)
#' ##lssplot_out<-plotlssVarReg(lssmodel, xlab="Time in seconds", ylab="Acceleration")
#' @export
#'

plotlssVarReg<-function(x, knot.lines=FALSE, xlab="x", ylab="y"){
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old.par))
  graphics::par(mfrow=c(2,2))
  if (is.logical(knot.lines)==FALSE){
    stop("knot.lines is not logical")
  }
  if (x$modeltype!="LSS model"){
    stop("This is not output from lssVarReg")
  }
  calc<-lss_calc(x)
  calc<-calc[order(calc[,2]),]
  graphics::plot(calc[,2],calc$`Predicted mean`,type="l",lwd=3,  col="red",xlab=xlab, ylab=ylab, ylim=c(min(calc$`Predicted mean`), max(calc$`Predicted mean`)), main="Predicted mean")
  if (x$locationmodel=="semi" && knot.lines==TRUE){
      bmean<-splines::bs(x$data[,3], df=(x$degree+x$knots.l), degree=x$degree)
      k <- attr(bmean, "knots")
      for (i in 1:(length(k))){
        graphics::abline(v=k[i], lty=2)
      }
  }
  graphics::plot(calc[,2],calc$`Predicted variance`,type="l",lwd=3,  col="blue",xlab=xlab, ylab="Variance", ylim=c(min(calc$`Predicted variance`), max(calc$`Predicted variance`)), main="Predicted variance")
  if (x$scale2model=="semi" && knot.lines==TRUE){
      bvar<-splines::bs(x$data[,3], df=(x$degree+x$knots.sc), degree=x$degree)
      k <- attr(bvar, "knots")
      for (i in 1:(length(k))){
        graphics::abline(v=k[i], lty=2)
      }
    }
  graphics::plot(calc[,2],calc$`Predicted skewness`,type="l",lwd=3,  col="dark green",xlab=xlab, ylab="Skew", ylim=c(min(calc$`Predicted skewness`), max(calc$`Predicted skewness`)), main="Predicted skewness")
  if (x$shapemodel=="semi" && knot.lines==TRUE){
    bvar<-splines::bs(x$data[,3], df=(x$degree+x$knots.sh), degree=x$degree)
    k <- attr(bvar, "knots")
    for (i in 1:(length(k))){
      graphics::abline(v=k[i], lty=2)
    }
  }
  n<-length(x$data[,1])
  chiq<-stats::qchisq(c(1:n)/(n+1), df=1)
  graphics::plot(chiq, sort(calc$stand.res2), main="Residual plot", pch=19, ylab="Squared standardised residuals", xlab="Chi-Square(1) quantiles")
  graphics::abline(0,1, col="red")
  return(calc)
}

