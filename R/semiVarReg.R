#' Semi parametric mean and variance regression
#'
#' \code{semiVarReg} performs semiparametric mean and variance regression. At the moment, this is only designed for a single covariate that is fit in the mean and variance models.
#' @param y Vector containing outcome data. Must be no missing data.
#' @param x Vector containing the covariate data. Must be no missing data.
#' @param meanmodel Specify the mean model to be fit: "zero" =zero mean model, "constant" = constant mean, "linear" = linear term with x covariate, "semi" = semiparametric spline (specify with $knots.m$)
#' @param varmodel Specify the variance model to be fit: "constant" = constant variance model, "linear" = linear term with x covariate, "semi" = semiparametric spline (specify with $knots.v$)
#' @param knots.m A integer indicating the number of internal knots to be fit in the mean model. Default is '2'. (Note that the knots are placed equidistantly over x.)
#' @param knots.v A integer indicating the number of internal knots to be fit in the variance model. Default is '2'. (Note that the knots are placed equidistantly over x.)
#' @param degree The degree of the splines fit in the mean and variance. Default is '2'.
#' @param mono.var Indicates whether the variance model is monotonic (only applicable for 'linear' or 'semi' variance models). Default is 'none' (no monotonic constraints). Options are 'inc' for increasing or 'dec' for decreasing.
#' @param maxit Number of maximum iterations for the EM algorithm, default 1000.
#' @param eps Very small number for the convergence criteria, default 1 times 10 power6
#' @return $semiVarReg$ returns an object of class "VarReg" which inherits some components from the class "glm".
#'
#' A list of the results from the algorithm, including conv, reldiff, information criterion and mean and variance estimates.
#'@examples
#'data(lidar)
#'linmodel<-semiVarReg(lidar$logratio, lidar$range, meanmodel="linear", varmodel="linear")
#'semimodel<-semiVarReg(lidar$logratio, lidar$range, meanmodel="semi", varmodel="semi", knots.m=4, knots.v=2)
#semimodel_inc<-semiVarReg(lidar$logratio, lidar$range, meanmodel="semi", varmodel="semi", knots.m=4, knots.var=2, mono)
#'@export


semiVarReg<-function(y,x, meanmodel=c("zero", "constant", "linear", "semi"), varmodel=c("constant", "linear", "semi"), knots.m=2, knots.v=2, degree=2, mono.var=c("none", "inc", "dec"), eps=1e-6, maxit=1000){
  meanmodel<-match.arg(meanmodel)
  varmodel<-match.arg(varmodel)
  mono.var<-match.arg(mono.var)
  if (length(y)!=length(x)){
    stop("Error: x and y not same length")
  }
  if ((meanmodel=="semi" || varmodel=="semi") && (length(knots.m)>1 || length(knots.v)>1)){
    stop("Error: knots.m or knots.v greater than length 1")
  }
  para.space<-"all"
  alldat<-data.frame(y,x)
  colnames(alldat)<-c(deparse(substitute(y)), deparse(substitute(x)))
  if (meanmodel=="semi"){
    bmean<-bs(x, df=(degree+knots.m), degree=degree)
    colnames(bmean) <- paste(paste("M_Knt",knots.m, sep = ""), paste("Base", colnames(bmean), sep=""), sep = "_")
    alldat<-data.frame(alldat, bmean)
    mean.ind<-which(colnames(alldat)%in%colnames(bmean))
  }else if (meanmodel=="linear"){
    mean.ind<-2
  }else if (meanmodel=="constant"){
    mean.ind<-0
  }else if (meanmodel=="zero"){
    mean.ind<-NULL
  }

  if (varmodel=="semi"){
  bvar<-bs(x, df=(degree+knots.v), degree=degree)
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
    var.ind<-2
    if (mono.var=="inc"){
      para.space<-"positive"
    }else if (mono.var=="inc"){
      para.space<-"negative"
    }
  }else if (varmodel=="constant"){
    var.ind<-FALSE
  }

  result<-linVarReg(alldat, var.ind=var.ind, mean.ind=mean.ind, para.space=para.space, eps=eps, maxit=maxit)
  model<-list(knots.m=knots.m, knots.v=knots.v, degree=degree, meanmodel=meanmodel, varmodel=varmodel)
  result<-c(model, result)
  return(result)
}