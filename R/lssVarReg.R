#' Semi parametric location, shape and scale regression
#'
#' \code{lssVarReg} performs a semiparametric location (\eqn{\xi} or xi), shape (\eqn{\nu} or nu) and scale (\eqn{\omega} or omega) regression model. Currently, this is only designed for a single covariate that is fit in the location, scale and shape models.
#' @param y Vector containing outcome data. Must be no missing data.
#' @param x Vector containing the covariate data, same length as \code{y}. Must be no missing data.
#' @param locationmodel Text to specify the location model to be fit. Options: \code{"constant"} = constant model (intercept only), \code{"linear"} = linear term with x covariate, \code{"semi"} = semiparametric spline (specify with \code{knots.l}).
#' @param scale2model Text to specify the scale^2 model to be fit. Options: \code{"constant"} = constant term only, \code{"linear"} = linear term with \code{x} covariate, \code{"semi"} = semiparametric spline (specify with \code{knots.sc})
#' @param shapemodel Text to specify the shape model to be fit. Options: \code{"constant"} = constant shape model, \code{"linear"} = linear term with x covariate, \code{"semi"} = semiparametric spline (specify with \code{knots.sh}).
#' @param knots.l Integer indicating the number of internal knots to be fit in the location model. Default is '2'. (Note that the knots are placed equidistantly over x.)
#' @param knots.sc Integer indicating the number of internal knots to be fit in the scale^2 model. Default is '2'. (Note that the knots are placed equidistantly over x.)
#' @param knots.sh Integer indicating the number of internal knots to be fit in the shape model. Default is '2'. (Note that the knots are placed equidistantly over x.)
#' @param degree Integer to indicate the degree of the splines fit in the location and scale.  Default is '2'.
#' @param mono.scale Text to indicate whether the scale2 model is monotonic. Default is \code{"none"} (no monotonic constraints). Options are \code{"inc"} for increasing or \code{"dec"} for decreasing. If this is chosen, the appropriate \code{para.space} is set autopmatically (\code{"positive"} for \code{inc}, \code{"negative"} for \code{dec}).
#' @param para.space Text to indicate the parameter space to search for scale2 parameter estimates. \code{"positive"} means only search positive parameter space, \code{"negative"} means search only negative parameter space and \code{"all"} means search all parameter spaces. Default is \code{all}.
#' @param location.init Vector of initial parameter estimates for the location model. Defaults to vector of 1's of appropriate length.
#' @param scale2.init Vector of initial parameter estimates for the scale^2 model. Defaults to vector of 1's of appropriate length.
#' @param shape.init Vector of initial parameter estimates for the shape model. Defaults to vector of 1's of appropriate length.
#' @param int.maxit Integer of maximum iterations for the internal location and scale EM algorithm. Default is 1000 iterations.
#' @param print.it Logical for printing progress of estimates through each iteration. Default is \code{FALSE}.
#' @param control List of control parameters for the algorithm. See \code{\link{VarReg.control}}.
#' @param ... arguments to be used to form the default control argument if it is not supplied
#' directly
#' @return
#' \code{lssVarReg} returns an object of class \code{"lssVarReg"}, which inherits most from class
#' \code{"VarReg"}. This object of class \code{lssVarReg} is a list of the following components:
#' \itemize{
#' \item \code{modeltype}: Text indicating the model that was fit, always "LSS model".
#' \item \code{locationmodel}, \code{scale2model}, \code{shapemodel}, \code{knots.l}, \code{knots.sc},
#' \code{knots.sh}, \code{degree},\code{mono.scale} : Returning the input variables as described above
#' \item\code{converged}: Logical argument indicating if convergence occurred.
#' \item\code{iterations}: Total iterations performed of the main algorithm (not including the
#'  internal EM algorithm).
#'  \item\code{reldiff}: the positive convergence tolerance that occured at the final iteration.
#'  \item\code{loglik}: Numeric variable of the maximised log-likelihood.
#'  \item\code{aic.c}: Akaike information criterion corrected for small samples
#'  \item\code{aic}: Akaike information criterion
#'  \item\code{bic}: Bayesian information criterion
#'  \item\code{hqc}: Hannan-Quinn information criterion
#'  \item\code{location}: Vector of the maximum likelihood estimates of the location parameters.
#'  \item\code{scale2}: Vector of the maximum likelihood estimates of the scale (squared) parameters.
#'  \item\code{shape}: Vector of the maximum likelihood estimates of the shape parameters.
#'  \item\code{data}: Dataframe containing the variables included in the model.
#'  }
#'
#'@seealso
#'  \code{\link{VarReg.control}}  \code{\link{plotlssVarReg}}
#'
#'@examples
# 'data(mcycle)
#' ## run a model with linear mean, linear variance and constant shape (not run):
#' ## lssmodel<-lssVarReg(mcycle$accel, mcycle$times,  locationmodel="linear", scale2model="linear",
#' ## shapemodel="constant",  maxit=10000)
#' @export


lssVarReg<-function(y, x, locationmodel=c("constant", "linear", "semi"),
                    scale2model=c("constant", "linear", "semi"),
                    shapemodel=c("constant", "linear"),
                    knots.l=2, knots.sc=2,knots.sh=2, degree=2,
                    mono.scale=c("none", "inc", "dec"),
                    para.space=c("all", "positive", "negative"),
                    location.init=NULL, scale2.init=NULL,shape.init=NULL, int.maxit=1000, print.it=FALSE, control=list(...), ...) {
  locationmodel<-match.arg(locationmodel)
  scale2model<-match.arg(scale2model)
  shapemodel<-match.arg(shapemodel)
  mono.scale<-match.arg(mono.scale)
  para.space<-match.arg(para.space)
  control<-do.call(VarReg.control, control)
  if (length(y)!=length(x)){
    stop("Error: x and y not same length")
  }
  if ((locationmodel=="semi" || scale2model=="semi") && (length(knots.l)>1 || length(knots.sc)>1)){
    stop("Error: knots.l or knots.sc greater than length 1")
  }
  n<-length(y)

  alldat<-data.frame(y, rep(1,n), x)
  colnames(alldat)<-c(deparse(substitute(y)),"mean.int", deparse(substitute(x)))
  if (locationmodel=="constant"){
    mean.ind<-c(2)
    if (is.null(location.init)==TRUE){
      location.init<-c(1)
    } else if (length(location.init)!=1){
      stop("Error: length of location.init doesnt match model requirements. (length of vector not correct for # of parameters)")
    }
  }else if (locationmodel=="linear"){
    mean.ind<-c(2,3)
    if (is.null(location.init)==TRUE){
      location.init<-c(1,1)
    } else if (length(location.init)!=2){
      stop("Error: length of location.init doesnt match model requirements. (length of vector not correct for # of parameters)")
    }
  }else if (locationmodel=="semi"){
    bmean<-splines::bs(x, df=(degree+knots.l), degree=degree)
    colnames(bmean) <- paste(paste("L_Knt",knots.l, sep = ""), paste("Base", colnames(bmean), sep=""), sep = "_")
    alldat<-data.frame(alldat, bmean)
    mean.ind<-c(2, which(colnames(alldat)%in%colnames(bmean)))
    if  (is.null(location.init)==TRUE){
      location.init<-rep(1, knots.l+3)
    } else if (length(location.init)!=(knots.l+3)){
      stop("Error: length of location.init doesnt match model requirements (length of vector not correct for # of parameters).")
    }
  }
  if (scale2model=="linear"){
    knots.sc<-NULL
    xvar<-as.data.frame(x)
    colnames(xvar)<-"Scale2_covariate"
    alldat<-data.frame(alldat, xvar)
    var.ind<-c(which(colnames(alldat)%in%colnames(xvar)))
    if (is.null(scale2.init)==TRUE){
      scale2.init<-c(1,1)
    } else if (length(scale2.init)!=2){
      stop("Error: length of scale2.init doesnt match model requirements (length of vector not correct for # of parameters).")
    }
  }else if (scale2model=="constant"){
    knots.sc<-NULL
    var.ind<-FALSE
    if (is.null(scale2.init)==TRUE){
    scale2.init<-1
    }else if (length(scale2.init)!=1){
      stop("Error: length of scale2.init doesnt match model requirements (length of vector not correct for # of parameters).")
    }
    }else if (scale2model=="semi"){
    bvar<-splines::bs(x, df=(degree+knots.sc), degree=degree)
    colnames(bvar) <- paste(paste("Sc_Knt",knots.sc, sep = ""), paste("Base", colnames(bvar), sep=""), sep = "_")
    if (mono.scale=="inc"){
      para.space="positive"
      bvar<-t(apply(bvar, 1, cumsum))
      colnames(bvar) <- paste(paste("INC_", colnames(bvar), sep=""), sep = "_")
    }else if (mono.scale=="dec"){
      para.space="negative"
      bvar<-t(apply(bvar[,rep(ncol(bvar):1)], 1, cumsum))
      bvar<-bvar[,rep(ncol(bvar):1)]
      colnames(bvar) <- paste(paste("DEC_", colnames(bvar), sep=""), sep = "_")
    }
    alldat<-data.frame(alldat, bvar)
    var.ind<-which(colnames(alldat)%in%colnames(bvar))
    if  (is.null(scale2.init)==TRUE){
      scale2.init<-rep(1, knots.sc+3)
    } else if (length(scale2.init)!=(knots.sc+3)){
      stop("Error: length of scale2.init doesnt match model requirements. (length of vector not correct for # of parameters)")
    }
  }
  if (shapemodel=="constant"){
    nuold<-0
    knots.sh<-NULL
    nu.ind<-NULL
    if (is.null(shape.init)==TRUE){
      shape.init=1
    }else if (length(shape.init)!=1){
      stop("Error: length of shape.init doesnt match model requirements (length of vector not correct for # of parameters).")
    }
  } else if (shapemodel=="linear"){
    nu.ind<-3
    if (is.null(shape.init)==TRUE){
      shape.init<-c(1,1)
    } else if (length(shape.init)!=2){
      stop("Error: length of shape.init doesnt match model requirements (length of vector not correct for # of parameters).")
    }
  }else if (shapemodel=="semi"){
    bshape<-splines::bs(x, df=(degree+knots.sh), degree=degree)
    colnames(bshape) <- paste(paste("Sh_Knt",knots.sh, sep = ""), paste("Base", colnames(bshape), sep=""), sep = "_")
    alldat<-data.frame(alldat, bshape)
    nu.ind<-which(colnames(alldat)%in%colnames(bshape))
    if  (is.null(shape.init)==TRUE){
      shape.init<-rep(1, knots.sh+3)
    } else if (length(shape.init)!=(knots.sh+3)){
      stop("Error: length of shape.init doesnt match model requirements. (length of vector not correct for # of parameters)")
    }
  }
  xiold<-location.init
  omega2old<-scale2.init
  nuold<-shape.init
    l<-loop_lss(alldat,xiold,omega2old,nuold,mean.ind, var.ind, nu.ind, para.space,maxit=control$maxit, eps=control$epsilon,int.maxit, print.it)
    mean<-l$fitted.xi
    variance<-unname(colSums(t(cbind(rep(1,n),alldat[1:n,var.ind]))*l$omega2new))
    nu<-unname(colSums(t(cbind(rep(1,n),alldat[1:n,nu.ind]))*l$nunew))
    d<-vector()
    for (i in 1:n){
      d[i]<-(sn::dsn(y[i], xi=l$fitted.xi[i], omega=sqrt(variance[i]), alpha=nu[i], log=TRUE))
    }
    loglik<-sum(d)
    param<-length(l$xinew)+length(l$omega2new)+length(l$nunew)
    ic<-criterion(n, loglik, param)

  out<-list(modeltype="LSS model", locationmodel=locationmodel, knots.l=knots.l, scale2model=scale2model, knots.sc=knots.sc, shapemodel=shapemodel, knots.sh=knots.sh, degree=degree, converged=l$conv, iterations=l$it,reldiff=l$reldiff, loglik=loglik, aic.c=ic$aicc, aic=ic$aic,bic=ic$bic, mono.scale=mono.scale, hqc=ic$hqc, location=l$xinew,scale2=l$omega2new,shape=l$nunew, data=alldat)
  class(out) <- c("lssVarReg")
  return(out)
}

