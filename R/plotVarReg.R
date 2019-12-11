#' Plots graphics for a mean and variance regression model
#'
#' \code{plotVarReg} to produce graphics for models fit in this package.
#' @param x Object of class \code{VarReg} (see \code{\link{semiVarReg}}).
#' @param knot.lines  Logical to indicate if knot lines should be shown on graphics
#' (if model is type "semi"). Default is \code{FALSE}
#' @param ci Logical indicate if 95\% CI should be shown on the plots. Default is \code{FALSE}
#' and \code{ci.type="im"}.
#' @param ci.type Text to indicate the type of CI to plot. Either \code{"im"} (information matrix) or \code{"boot"} (bootstrapped). Default is \code{"im"}.
#' @param bootreps Integer to indicate the number of bootstrap replications to be performed if \code{ci.type="boot"}. Default is \code{1000}.
#' @param xlab Text for the label to be placed on the \code{x} axis of graphics (covariate)
#' @param ylab Text for the label to be placed on the \code{y} axis of graphics (outcome)
#' @param control list of control parameters to be used in bootstrapping.
#' See \code{\link{VarReg.control}}.
#' @param ... arguments to be used to form the default control argument if it is not supplied
#' directly
#' @return This function returns a 2x2 plot, with slightly different plots given, depending on the outcome data. For uncensored data, the plots are:
#' \itemize{
#' \item the mean function over the \code{x}-variable, with or without 95\% CI,  and with or
#' without the knot lines indicated
#' \item the variance function over the \code{x}-variable, with or without 95\% CI and with or
#' without the knot lines indicated
#' \item a Q-Q plot of the residuals from the model
#' \item a histogram of the residuals from the model
#' }
#'  If the outcome data is censored, the last two plots are no longer appropriate.
#'  Given the censored residuals from the model, we can compare the squared standardised residuals
#'  (given in black) with their censoring indicator to the chi-squared distribution with one
#'  degree of freedom (given in red). This is one of the plots given for censored data, and the
#'  other is a plot of the data, coloured by the censoring status. The triangles with the point at
#'  the top are bottom censored and the triangles with the point at the bottom are top censored.
#'
#'@seealso \code{\link{semiVarReg}}, \code{\link{VarReg.control}}
#'
#'@examples
#'data(mcycle)
#'linmodel<-semiVarReg(mcycle$accel, mcycle$times, meanmodel="linear", varmodel="linear",
#'maxit=10000)
#'plotVarReg(linmodel)
#'plotVarReg(linmodel, ci=TRUE, ci.type="im", ylab="Range", xlab="Time in seconds")
#'##not run
#'##plotVarReg(linmodel, ci=TRUE, ci.type="boot", bootreps=10,ylab="Acceleration",
#'##xlab="Time in seconds")
#'
#'##not run
#'##semimodel<-semiVarReg(mcycle$accel, mcycle$times, meanmodel="semi", varmodel="semi",
#'##knots.m=4, knots.v=2, maxit=10000)
#'##plotVarReg(semimodel, ci=TRUE, ci.type="boot",bootreps=10,ylab="Acceleration",
#'##xlab="Time in seconds", maxit=10000)
#' @export


plotVarReg<-function(x, knot.lines=FALSE, ci=FALSE, ci.type=c("im", "boot"),bootreps=1000, xlab="x", ylab="y", control=list(...), ...){
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old.par))
  graphics::par(mfrow=c(2,2))
  if (is.logical(knot.lines)==FALSE){
    stop("Error: knot.lines is not logical")
  }
  if (is.logical(ci)==FALSE){
    stop("Error: ci is not logical")
  }
  ci.type<-match.arg(ci.type)
if (x$boundary==TRUE && ci.type=="im"){
  print("Warning: ci.type changed to 'boot' as boundary==TRUE")
  ci.type<-"boot"
}
  long<-seq(min(x$data[,2]), max(x$data[,2]), by=((max(x$data[,2])-min(x$data[,2]))/999))
  if (ci==FALSE){
    se_noCI<-seVarReg(x, boot=FALSE, vector.mean = long, vector.variance = long)
    graphics::plot(se_noCI$mean.output[,1],se_noCI$mean.output[,2],type="l",lwd=3,  col="red",xlab=xlab, ylab=ylab, ylim=c(min(se_noCI$mean.output[,2]), max(se_noCI$mean.output[,2])), main="Predicted mean")
    if (x$meanmodel=="semi" && knot.lines==TRUE){
      bmean<-splines::bs(x$data[,2], df=(x$degree+x$knots.m), degree=x$degree)
      k <- attr(bmean, "knots")
      for (i in 1:(length(k))){
        graphics::abline(v=k[i], lty=2)
      }
    }
    graphics::plot(se_noCI$variance.outputs[,1],se_noCI$variance.outputs[,2],type="l",lwd=3,  col="blue",xlab=xlab, ylab="Variance", ylim=c(min(se_noCI$variance.outputs[,2]), max(se_noCI$variance.outputs[,2])), main="Predicted variance")
    if (x$varmodel=="semi" && knot.lines==TRUE){
      bvar<-splines::bs(x$data[,2], df=(x$degree+x$knots.v), degree=x$degree)
      k <- attr(bvar, "knots")
      for (i in 1:(length(k))){
        graphics::abline(v=k[i], lty=2)
      }
    }
  } else if (ci==TRUE && ci.type=="im"){
    print("CI=true, type=information matrix")
    se_CI<-seVarReg(x, boot=FALSE, vector.mean = long, vector.variance = long)
    if (x$meanmodel=="semi"){
      print("Meanmodel='semi' so 95% CI cannot be given by information matrix")
      graphics::plot(se_CI$mean.output[,1],se_CI$mean.output[,2],type="l",lwd=3,  col="red",xlab=xlab, ylab=ylab, ylim=c(min(se_CI$mean.output[,2]), max(se_CI$mean.output[,2])), main="Predicted mean")
      if (knot.lines==TRUE){
        bmean<-splines::bs(x$data[,2], df=(x$degree+x$knots.m), degree=x$degree)
        k <- attr(bmean, "knots")
        for (i in 1:(length(k))){
          graphics::abline(v=k[i], lty=2)
        }
      }
    }else if (x$meanmodel=="zero"){
      print("Meanmodel='zero' so 95% CI cannot be given by information matrix")
      graphics::plot(se_CI$mean.output[,1],se_CI$mean.output[,2],type="l",lwd=3,  col="red",xlab=xlab, ylab=ylab, ylim=c(min(se_CI$mean.output[,2]), max(se_CI$mean.output[,2])), main="Predicted mean")
    }else {
      mean.lci<-se_CI$mean.outputs[,2]-1.96*se_CI$mean.outputs[,3]
      mean.uci<-se_CI$mean.outputs[,2]+1.96*se_CI$mean.outputs[,3]

      graphics::plot(se_CI$mean.outputs[,1],se_CI$mean.outputs[,2],type="n",xlab=xlab, ylab=ylab,
           ylim=c(min(mean.lci), max(mean.uci)), main="Predicted mean")
      graphics::polygon(c(se_CI$mean.outputs[,1], rev(se_CI$mean.outputs[,1])), c(mean.lci,rev(mean.uci)), col = "lightgray", border = NA)
      graphics::lines(se_CI$mean.outputs[,1],se_CI$mean.outputs[,2], lwd=3,  col="red")
    }
    if (x$varmodel=="semi"){
      print("Varmodel='semi' so 95% CI cannot be given by information matrix")
      graphics::plot(se_CI$variance.outputs[,1],se_CI$variance.outputs[,2],type="l",lwd=3,  col="blue",xlab=xlab, ylab="Variance", ylim=c(min(se_CI$variance.outputs[,2]), max(se_CI$variance.outputs[,2])), main="Predicted variance")
    }else{
      var.lci<-se_CI$variance.outputs[,2]-1.96*se_CI$variance.outputs[,3]
      var.uci<-se_CI$variance.outputs[,2]+1.96*se_CI$variance.outputs[,3]

      graphics::plot(se_CI$variance.outputs[,1],se_CI$variance.outputs[,2],type="n",xlab=xlab, ylab="Variance", ylim=c(min(var.lci), max(var.uci)), main="Predicted variance")
      graphics::polygon(c(se_CI$variance.outputs[,1], rev(se_CI$variance.outputs[,1])), c(var.lci,rev(var.uci)), col = "lightgray", border = NA)
      graphics::lines(se_CI$variance.outputs[,1],se_CI$variance.outputs[,2], lwd=3,  col="blue")
    }


  }else if (ci==TRUE && ci.type=="boot"){
    print("CI=true, type=bootstrapped")
    se_boot<-seVarReg(x, boot=TRUE,bootreps=bootreps,control=control, vector.mean = long, vector.variance = long)
    if (x$meanmodel!="zero"){
      mean.lci<-se_boot$mean.outputs[,5]
      mean.uci<-se_boot$mean.outputs[,6]
      graphics::plot(se_boot$mean.outputs[,1],se_boot$mean.outputs[,2],type="n",xlab=xlab, ylab=ylab,
           ylim=c(min(mean.lci), max(mean.uci)), main="Predicted mean")
      graphics::polygon(c(se_boot$mean.outputs[,1], rev(se_boot$mean.outputs[,1])), c(mean.lci,rev(mean.uci)), col = "lightgray", border = NA)
      graphics::lines(se_boot$mean.outputs[,1],se_boot$mean.outputs[,2], lwd=3,  col="red")
      if (x$meanmodel=="semi" && knot.lines==TRUE){
        bmean<-splines::bs(x$data[,2], df=(x$degree+x$knots.m), degree=x$degree)
        k <- attr(bmean, "knots")
        for (i in 1:(length(k))){
          graphics::abline(v=k[i], lty=2)
          }
      }
    }else if (x$meanmodel=="zero"){
      print("No mean plot with 95% CI for zero mean model")
      graphics::plot(se_boot$mean.output[,1],se_boot$mean.output[,2],type="l",lwd=3,  col="red",xlab=xlab, ylab=ylab, ylim=c(min(se_boot$mean.output[,2]), max(se_boot$mean.output[,2])), main="Predicted mean")
    }

    var.lci<-se_boot$variance.outputs[,5]
    var.uci<-se_boot$variance.outputs[,6]
    graphics::plot(se_boot$variance.outputs[,1],se_boot$variance.outputs[,2],type="n",xlab=xlab, ylab="Variance", ylim=c(min(var.lci), max(var.uci)), main="Predicted variance")
    graphics::polygon(c(se_boot$variance.outputs[,1], rev(se_boot$variance.outputs[,1])), c(var.lci,rev(var.uci)), col = "lightgray", border = NA)
    graphics::lines(se_boot$variance.outputs[,1],se_boot$variance.outputs[,2], lwd=3,  col="blue")
    if (x$meanmodel=="semi" && knot.lines==TRUE){
      bmean<-splines::bs(x$data[,2], df=(x$degree+x$knots.m), degree=x$degree)
      k <- attr(bmean, "knots")
      for (i in 1:(length(k))){
        graphics::abline(v=k[i], lty=2)
      }
    }
  }
  se_res<-seVarReg(x, boot=FALSE)
  stand.res<-(x$data[,1]-se_res$mean.outputs[,2])/sqrt(se_res$variance.outputs[,2])
  if (is.null(x$cens.ind)==TRUE){
    stats::qqnorm(stand.res, pch=19)
    stats::qqline(stand.res, col=2)
    graphics::hist(stand.res,  main="Histogram of residuals", xlab="Standardised residuals")
  }else {
    censor<-vector()
    censor[x$data[,x$cens.ind]==0]<-1
    censor[x$data[,x$cens.ind]==-1]<-0
    censor[x$data[,x$cens.ind]==1]<-0
    censor2<-x$data[,x$cens.ind]
    stand.res2<-stand.res**2
    x1 <- seq(0,max(stand.res2),by=.5)
    y <- stats::pchisq(x1,df=1)
    graphics::plot(x1,y, col="red", type="l", lwd=1, xlab="", ylab="", main = "Censored residuals")
    graphics::lines(survival::survfit(survival::Surv(stand.res2, censor==1)~1), col='black',lwd=2,  fun="event", conf.int=FALSE)
    graphics::plot(x$data[,2],x$data[,1],main="Censored obs in red", ylab=ylab, xlab=xlab,pch=c(2,19,6)[censor2+2],   col=c('red','black','red')[censor2+2] )
    }
}

