#' Plots graphics for a mean and variance regression model
#'
#' \code{plotVarReg} to produce graphics for models fit in this package.
#' @param x Object of class VarReg.
#' @param knot.lines  TRUE/FALSE (Default=TRUE) to show the knot lines on the graphics.
#' @param CI TRUE/FALSE (Default=FALSE) to show the 95\% CI on the plots.
#' @return returns a series of graphics.
#'
#'@examples
#'data(lidar)
#'linmodel<-semiVarReg(lidar$logratio, lidar$range, meanmodel="linear", varmodel="linear")
#'test<-plotVarReg(linmodel)
#' @export


plotVarReg<-function(x, knot.lines=FALSE, CI=FALSE, ci.type=c("im", "boot"), xlab="x", ylab="y"){
  if (is.logical(knot.lines)==FALSE){
    stop("Error: knot.lines is not logical")
  }
  if (is.logical(CI)==FALSE){
    stop("Error: CI is not logical")
  }
  ci.type<-match.arg(ci.type)

  n<-length(x$data[,1])
  if (x$meanmodel=="zero"){
    pred.mean<-NULL
    CI=FALSE
  }else if (x$meanmodel=="constant" ||x$meanmodel=="linear"){
    d<-cbind(rep(1,n), x$data[,x$mean.ind])
    pred.mean<-colSums(t(d)*x$mean)
    newx<-seq(min(x$data[,x$mean.ind]), max(x$data[,x$mean.ind]), by =(max(x$data[,x$mean.ind])- min(x$data[,x$mean.ind]))/999)
    grid.d<-cbind(rep(1,1000), newx)
    grid.mean<-colSums(t(grid.d)*x$mean)
  }else if (x$meanmodel=="semi"){
    d<-cbind(rep(1,n), x$data[,x$mean.ind])
    pred.mean<-colSums(t(d)*x$mean)
    knots.mean<-splines::bs(x$data[,2], df=(x$degree+x$knots.m), degree=x$degree)
    k <- attr(knots.mean, "knots")
    b <- attr(knots.mean, "Boundary.knots")
    newx<-seq(b[1], b[2], by=(b[2]-b[1])/999)
    knots2<-splines::bs(newx, df=(length(k)+2), degree=x$degree, knots=k, Boundary.knots=b)
    grid.d<-cbind(rep(1,1000), knots2)
    grid.mean<-colSums(t(grid.d)*x$mean)
  }

  if (CI==FALSE){
    plot(newx,grid.mean,type="l",lwd=3,  col="red",xlab=xlab, ylab=ylab, ylim=c(min(grid.mean), max(grid.mean)), main="Predicted mean")
  } else if (CI==TRUE && ci.type=="im"){
    se<-seVarReg(x, boot=FALSE)
    mean.var<-vector(length=length(pred.mean))

    for (i in 1:ncol(grid.d)){
      for (j in 1:ncol(grid.d)){
        mean.var<-mean.var+grid.d[,i]*grid.d[,j]*se$mean.im[i,j]
        grid.mean.var<-grid.mean.var+grid.d[,i]*grid.d[,j]*se$mean.im[i,j]
      }
    }
    mean.se<-
    mean.lci<-pred.mean-1.96*sqrt(mean.var)
    mean.uci<-pred.mean+1.96*sqrt(mean.var)
  }
  if (meanmodel=="semi" && knot.lines==TRUE){
    for (i in 1:(length(k))){
      abline(v=k[i], lty=2)
    }
  }
}
#
#   par(mfrow=c(1,1))
#   ##create grid of mean and variance of mean
#   if (is.null(knots.mean)==FALSE){
#     k <- attr(knots.mean, "knots")
#     b <- attr(knots.mean, "Boundary.knots")
#     x<-seq(b[1], b[2], by=(b[2]-b[1])/1000)
#     knots2<-bs(x, df=(length(k)+2), degree=2, knots=k, Boundary.knots=b)
#     knots.m<-data.frame(rep(1, length(x)), knots2)
#   }else if (length(outmodel$mean.ind)==1){
#     k=NULL
#     x<-seq(min(outmodel$dat[outmodel$mean.ind,]), max(outmodel$dat[outmodel$mean.ind,]), by =(max(outmodel$dat[outmodel$mean.ind,])- min(outmodel$dat[outmodel$mean.ind,])/1000))
#     knots.m<-data.frame(rep(1, length(x)), x)
#   }
#
#
#
#   ##pred mean
#   pred.mean<-colSums(t(knots.m)*outmodel$mean[,1])
#
#   ##pred var of mean
#   mean.var<-vector(length=length(pred.mean))
#
#   for (i in 1:ncol(knots.m)){
#     for (j in 1:ncol(knots.m)){
#       mean.var<-mean.var+knots.m[,i]*knots.m[,j]*outmodel$im.beta[i,j]
#     }
#   }
#   mean.lci<-pred.mean-1.96*sqrt(mean.var)
#   mean.uci<-pred.mean+1.96*sqrt(mean.var)
#
#   ##plot mean
#   plot(x,pred.mean,type="l",lwd=3,  col="red",xlab=xlab, ylab=ylab, ylim=c(min(mean.lci), max(mean.uci)), main="Predicted mean")
#   lines(x,mean.lci, col="black", lwd=3)
#   lines(x,mean.uci, col="black", lwd=3)
#   if (lines==TRUE && is.null(k)!=TRUE){
#     for (i in 1:(length(k))){
#       abline(v=k[i])
#     }
#   }
#
#   ##grid for pred var
#   if (is.null(knots.var)==FALSE){
#     k.v <- attr(knots.var, "knots")
#     b <- attr(knots.var, "Boundary.knots")
#     x<-seq(b[1], b[2], by=(b[2]-b[1])/1000)
#
#     knots2<-bs(x, df=(length(k.v)+2), degree=2, knots=k.v, Boundary.knots=b)
#     knots.v<-data.frame(rep(1, length(x)), knots2)
#   }else if (length(outmodel$var.ind)==1){
#     k.v=NULL
#     x<-seq(min(outmodel$dat[outmodel$var.ind,]), max(outmodel$dat[outmodel$var.ind,]), by =(max(outmodel$dat[outmodel$var.ind,])- min(outmodel$dat[outmodel$var.ind,])/1000))
#     knots.v<-data.frame(rep(1, length(x)), x)
#   }
#
#   ##pred var
#   pred.var<-colSums(t(knots.v)*outmodel$variance[,1])
#   var.var<-vector(length=length(pred.var))
#
#
#   for (i in 1:ncol(knots.v)){
#     for (j in 1:ncol(knots.v)){
#       var.var<-var.var+knots.v[,i]*knots.v[,j]*outmodel$im.alpha[i,j]
#     }
#   }
#   var.lci<-pred.var-1.96*sqrt(var.var)
#   var.uci<-pred.var+1.96*sqrt(var.var)
#
#   ##plot predvar
#   plot(x,pred.var, col="darkgreen", type="l",lwd=3, ylim=c(min(var.lci), max(var.uci)), xlab=xlab, ylab="Variance", main="Predicted Variance")
#   lines(x,var.lci, col="black", lwd=3)
#   lines(x,var.uci, col="black", lwd=3)
#   if (lines==TRUE && is.null(k.v)!=TRUE){
#     for (i in 1:(length(k.v))){
#       abline(v=k.v[i])
#     }
#   }
#   ##plot SD
#   sd.lci<-sqrt(pred.var)-1.96*(1/sqrt(4*pred.var))*sqrt(var.var)
#   sd.uci<-sqrt(pred.var)+1.96*(1/sqrt(4*pred.var))*sqrt(var.var)
#
#
#   plot(x,sqrt(pred.var), col="darkgreen", type="l",lwd=3,  ylim=c(min(sd.lci), max(sd.uci)), xlab=xlab, ylab="Stddev", main="Predicted Stddev")
#   lines(x,sd.lci, col="black", lwd=3)
#   lines(x,sd.uci, col="black", lwd=3)
#   if (lines==TRUE && is.null(k.v)!=TRUE){
#     for (i in 1:(length(k.v))){
#       abline(v=k.v[i])
#     }
#   }
#   ##residuals
#   pred.mean.dat<-colSums(t(cbind(rep(1, length(outmodel$data[,1])),knots.mean))*outmodel$mean[,1])
#   pred.var.dat<-colSums(t(cbind(rep(1, length(outmodel$data[,1])),knots.var))*outmodel$variance[,1])
#   stand.res<-(outmodel$data[,1]-pred.mean.dat)/sqrt(pred.var.dat)
#   par(mfrow=c(1,2))
#   qqnorm(stand.res, pch=19)
#   qqline(stand.res, col=2)
#   plot(density(stand.res),  main="Histogram of residuals")


