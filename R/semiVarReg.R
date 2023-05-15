#' Semi parametric mean and variance regression
#'
#' \code{semiVarReg} performs semi-parametric mean and variance regression models. Currently, this is
#' only designed for a single covariate that is fit in the mean and variance models.
#'@param y Vector containing outcome data. Must be no missing data and any censored values must
#' be set to the limits of detection.
#'@param x Vector containing the covariate data. Must be no missing data and same length as \code{y}.
#'@param cens.ind Vector containing the censoring indicator, if applicable. There must be no missing
#' data contained in the vector and this vector should be the same length as \code{y}.
#' \code{"0"} values indicate uncensored data, \code{"1"} indicates right, or upper, censoring and
#' \code{"-1"} indicates left, or lower, censoring. The default is \code{NULL} which indicates there
#' is no censored data.
#'@param meanmodel Text to specify the mean model to be fit to the data. The possible inputs are
#'  \code{"zero"}, \code{"constant"}, \code{"linear"} or \code{"semi"}. \code{"semi"}
#'  indicates a semi-parametric spline model, with the number of internal knots specified in
#'  \code{knots.m}.
#'@param mean.intercept Logical argument  to indicate if the mean model is to include an intercept
#' term. This option is only available in the censored mean model, and the default=\code{TRUE}.
#'@param varmodel Text to specify the variance model to be fit to the data. The possible inputs are
#' \code{"constant"}, \code{"linear"} or \code{"semi"}. \code{"semi"} indicates a semi-parametric
#' B-spline model, with the number of internal knots specified in \code{knots.v}.
#'@param knots.m Integer indicating the number of internal knots to be fit in the semi-parametric
#' mean model. Knots are placed equidistantly over the covariate. The default value is \code{2}.
#'@param knots.v Integer indicating the number of internal knots to be fit in the semi-parametric
#'  variance model. Knots are placed equidistantly over the covariate. The default value is \code{2}.
#'@param degree Integer indicating the degree of the splines fit in the mean and the variance models.
#' The default value is \code{2}.
#' @param mono.var Text to indicate whether the variance model is monotonic. Note that this is not
#'  available for the \code{"constant"} variance model. Options are \code{"none"}, \code{"inc"} or
#' \code{"dec"}, with the default=\code{"none"}. \code{"Inc"} indicates increasing monotonic and
#' \code{"dec"} indicates decreasing monotonic. If the variance model is linear, the parameter
#' space is constrained (positive for increasing and negative for decreasing). For semi-parametric
#' variance models, the appropriate monotonic B-splines are fit in the semi-parametric variance model.
#'@param para.space Text to indicate the parameter space to search for scale2 parameter estimates.
#' \code{"positive"} means only search positive parameter space, \code{"negative"} means search only
#'  negative parameter space and \code{"all"} means search all parameter spaces. Default is \code{all}.
#' @param control list of control parameters. See \code{\link{VarReg.control}}.
#' @param ... arguments to be used to form the default control argument if it is not supplied
#' directly
#' @return \code{semiVarReg} returns an object of class \code{"VarReg"} which inherits some components from the class \code{"glm"}. This object of class \code{"VarReg"} is a list containing the following components:
#' \itemize{
#'  \item\code{modeltype}: Text indicating the model that was fit, indicating if a censored approach or an uncensored approach was performed.
#'  \item\code{knots.m}, \code{knots.v}, \code{degree}, \code{meanmodel}, \code{varmodel}: Returning the input variables as described above
#'  \item\code{converged}: Logical argument indicating if convergence occurred.
#'  \item\code{iterations}: Total iterations performed.
#'  \item\code{reldiff}: the positive convergence tolerance that occurred at the final iteration.
#'  \item\code{loglik}: Numeric variable of the maximised log-likelihood.
#'  \item\code{boundary}: Logical argument indicating if the MLE is on the boundary of the parameter space.
#'  \item\code{aic.c}: Akaike information criterion corrected for small samples
#'  \item\code{aic}: Akaike information criterion
#'  \item\code{bic}: Bayesian information criterion
#'  \item\code{hqc}: Hannan-Quinn information criterion
#'  \item\code{mean.ind}: Vector of integer(s) indicating the column number(s) in the dataframe
#'  \code{data} that were fit in the mean model.
#'  \item\code{mean}: Vector of the maximum likelihood estimates of the mean parameters.
#'  \item \code{var.ind}: Vector of integer(s) indicating the column(s) in the dataframe
#'  \code{data} that were fit in the variance model.
#'  \item\code{variance}: Vector of the maximum likelihood estimates of the variance parameters.
#'  \item\code{cens.ind}: Integer indicating the column in the dataframe \code{data} that
#'  corresponds to the censoring indicator.
#'  \item\code{data}: Dataframe containing the variables included in the model.}
#'
#'
#'@examples
#'data(mcycle)
#'## run a model with linear mean and linear variance:
#'linmodel<-semiVarReg(mcycle$accel, mcycle$times, meanmodel="linear", varmodel="linear",
#'  maxit=10000)
#'## run a model with semi-parametric mean (4 internal knots) and semi-parametric variance (2 knots):
#'##not run
#'##semimodel<-semiVarReg(mcycle$accel, mcycle$times, meanmodel="semi", varmodel="semi",
#'##knots.m=4, knots.v=2, maxit=10000)
#'## run a model with semi-parametric mean (4 internal knots) and semi-parametric monotonic
#'## variance (2 knots):
#'## not run
#'##semimodel_inc<-semiVarReg(mcycle$accel, mcycle$times, meanmodel="semi", varmodel="semi",
#'##knots.m=4, knots.v=2, mono.var="inc")
#'@export


semiVarReg<-function(y,x, cens.ind=NULL, meanmodel=c("zero", "constant", "linear", "semi"), mean.intercept=TRUE, varmodel=c("constant", "linear", "semi"), knots.m=2, knots.v=2, degree=2, mono.var=c("none", "inc", "dec"),para.space=c("all", "positive", "negative"), control=list(...), ...){
  meanmodel<-match.arg(meanmodel)
  varmodel<-match.arg(varmodel)
  mono.var<-match.arg(mono.var)
  para.space<-match.arg(para.space)
  control<-do.call(VarReg.control, control)
  if (length(y)!=length(x)){
    stop("Error: x and y not same length")
  }
  if ((meanmodel=="semi" || varmodel=="semi") && (length(knots.m)>1 || length(knots.v)>1)){
    stop("Error: knots.m or knots.v greater than length 1")
  }
  alldat<-data.frame(y,x)
  colnames(alldat)<-c(deparse(substitute(y)), deparse(substitute(x)))
  if (meanmodel=="semi"){
    bmean<-splines::bs(x, df=(degree+knots.m), degree=degree)
    colnames(bmean) <- paste(paste("M_Knt",knots.m, sep = ""), paste("Base", colnames(bmean), sep=""), sep = "_")
    alldat<-data.frame(alldat, bmean)
    mean.ind<-which(colnames(alldat)%in%colnames(bmean))
  }else if (meanmodel=="linear"){
    mean.ind<-2
    knots.m<-NULL
  }else if (meanmodel=="constant"){
    mean.ind<-0
    knots.m<-NULL
  }else if (meanmodel=="zero"){
    mean.ind<-NULL
    knots.m<-NULL
  }

  if (varmodel=="semi"){
  bvar<-splines::bs(x, df=(degree+knots.v), degree=degree)
  colnames(bvar) <- paste(paste("V_Knt",knots.v, sep = ""), paste("Base", colnames(bvar), sep=""), sep = "_")
    if (mono.var=="inc"){
      para.space<-"positive"
      bvar<-t(apply(bvar, 1, cumsum))
      colnames(bvar) <- paste(paste("V_Knt",knots.v, sep = ""), paste("INC_Base", colnames(bvar), sep=""), sep = "_")
    }else if (mono.var=="dec"){
      bvar<-t(apply(bvar[,rep(ncol(bvar):1)], 1, cumsum))
      bvar<-bvar[,rep(ncol(bvar):1)]
      colnames(bvar) <- paste(paste("V_Knt",knots.v, sep = ""), paste("DEC_Base", colnames(bvar), sep=""), sep = "_")
      para.space<-"negative"
    }
  alldat<-data.frame(alldat, bvar)
  var.ind<-which(colnames(alldat)%in%colnames(bvar))
  }else if (varmodel=="linear"){
    knots.v<-NULL
    var.ind<-2
    if (mono.var=="inc"){
      para.space<-"positive"
    }else if (mono.var=="dec"){
      para.space<-"negative"
    }
  }else if (varmodel=="constant"){
    knots.v<-NULL
    var.ind<-FALSE
  }

  if (is.null(cens.ind)==TRUE){
    result<-linVarReg(alldat, var.ind=var.ind, mean.ind=mean.ind, para.space=para.space, control=control)
    type<-"Mean and Variance regression"
  }else {
    alldat<-cbind(alldat,cens.ind)
    c.ind<-which(colnames(alldat)%in%("cens.ind"))
    result<-censlinVarReg(dat=alldat, var.ind=var.ind, mean.ind=mean.ind,mean.intercept=mean.intercept, cens.ind=c.ind, para.space=para.space, control=control)
    type<-c("Censored Mean and Variance regression with mean.intercept=",mean.intercept)
  }

  model<-list(modeltype=type, knots.m=knots.m, knots.v=knots.v, degree=degree, meanmodel=meanmodel, varmodel=varmodel)
  result<-c(model, result)
  class(result) <- c("VarReg")
  return(result)
}