#' Calculations for SN
#'
#' \code{lss_calc} performs calculations for transforming SN data (location, scale and shape) to mean, variance and skew. This function is utilised by other, higher level functions.
#' @param x Object of class lssVarReg (output from \code{lssVarReg}).
#' @return dataframe containing:
#' \itemize{
#' \item\code{y}: y variable
#' \item\code{x}: x variable
#' \item\code{eta}: \eqn{\eta} or fitted location estimates
#' \item\code{omega}: \eqn{\omega} or fitted scale estimates
#' \item\code{shape}: \eqn{\alpha} or fitted shape estimates
#'  \item\code{predicted mean}: fitted mean estimates
#'  \item\code{predicted variance}: fitted variance estimates
#'  \item\code{Predicted skewness}: fitted skewness estimates
#'  \item\code{stand.res2}: Squared standardised residuals
#'}
#' @export

lss_calc<-function(x){
  n<-length(x$data[,1])
  if (x$locationmodel=="semi"){
    bmean<-splines::bs(x$data[,3], df=(x$degree+x$knots.l), degree=x$degree)
    eta<-colSums(t(cbind(rep(1,n),bmean))*x$location)
  }else if (x$locationmodel=="linear"){
    eta<-colSums(t(cbind(rep(1,n),x$data[,3]))*x$location)
  }
  if (x$scale2model=="linear"){
    omega<-sqrt(colSums(t(cbind(rep(1,n),x$data[,3]))*x$scale2))
  }else if (x$scale2model=="constant"){
    omega<-sqrt(colSums(t(cbind(rep(1,n)))*x$scale2))
  }else if (x$scale2model=="semi"){
    bvar<-splines::bs(x$data[,3], df=(x$degree+x$knots.sc), degree=x$degree)
    if (x$mono.scale=="inc"){
      bvar<-t(apply(bvar, 1, cumsum))
      }else if (x$mono.scale=="dec"){
      bvar<-t(apply(bvar[,rep(ncol(bvar):1)], 1, cumsum))
      bvar<-bvar[,rep(ncol(bvar):1)]
      }
    omega<-sqrt(colSums(t(cbind(rep(1,n),bvar))*x$scale2))
  }
  if (x$shapemodel=="constant"){
    shape<-rep(x$shape, n)
  }else if (x$shapemodel=="linear"){
    shape<-colSums(t(cbind(rep(1,n),x$data[,3]))*x$shape)
  }

  delta<-shape/(sqrt(1+shape**2))
  pred.mean<-eta+omega*delta*sqrt(2/pi)
  pred.var<-omega**2*(1-(2*delta**2)/pi)
  skewness<-((4-pi)/2)*((delta*sqrt(2/pi))**3/((1-(2*delta**2)/pi)**(3/2)))
  stand.res2<- (x$data[,1]-eta)**2/(omega)**2

  out<-data.frame(x$data[,1], x$data[,3], eta,omega,shape, pred.mean, pred.var, skewness, stand.res2)
  colnames(out)<-c(colnames(x$data)[1], colnames(x$data)[3], "eta", "omega", "shape", "Predicted mean", "Predicted variance", "Predicted skewness", "stand.res2")
  return(out)
}