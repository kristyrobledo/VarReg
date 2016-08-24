#' Searches for best semi parametric mean and variance regression model
#'
#' \code{searchVarReg} performs multiple semiparametric mean and variance regression models for a covariate of interest. The best model is chosen based on the information criterion of preference (\code{"selection"}). At the moment, this is only designed for a single covariate that is fit in both the mean and variance models.
#' @param y Vector containing outcome data. Must be no missing data.
#' @param x Vector containing the covariate data. Must be no missing data.
#' @param maxknots.m  A integer indicating the maximum number of internal knots to be fit in the mean model. Default is '3'. (Note that the knots are placed equidistantly over x.)
#' @param maxknots.v A integer indicating the maximum number of internal knots to be fit in the variance model. Default is '3'. (Note that the knots are placed equidistantly over x.)
#' @param degree The degree of the splines fit in the mean and variance. Default is '2'.
#' @param mono.var Indicates whether the variance model is monotonic (only applied to 'linear' or splines variance models). Default is 'none' (no monotonic constraints). Options are 'inc' for increasing or 'dec' for decreasing.
#' @param selection Indicates which information criteria is to be used for the selection of the best model. Choices are 'AIC', 'AICc', 'HQC', 'BIC'.
#' @param control  list of control parameters. See \code{\link{VarReg.control}}.
#' @return \code{searchVarReg} returns an list, where one of the components \code{best.model} is of class \code{"VarReg"} which inherits some components from the class \code{"glm"}. The list also contains the matrix of AIC, AICc, BIC and HQC from all of the models fit to the data.
#'@details A matrix of models are performed, of increasing complexity. Mean models start at a zero mean model, then constant mean, linear, 0 internal knots, etc, up to a maximum internal knots as specified in \code{maxknots.m}. Variance models start at constant variance, linear variance, 0 internal knots, etc, up to max internal knots as specified in \code{maxknots.v}.

#'@examples
#'data(lidar)
#'find<-searchVarReg(lidar$logratio, lidar$range, maxknots.v=3, maxknots.m=3, selection="HQC", maxit=100)

#'@export

searchVarReg<-function(y,x, maxknots.m=3, maxknots.v=3, degree=2, mono.var=c("none", "inc", "dec"), selection=c("AIC", "AICc", "HQC", "BIC"), control=list(...), ...){
  selection<-match.arg(selection)
  if (length(y)!=length(x)){
    stop("Error: x and y not same length")
  }
  if (length(maxknots.m)>1 || length(maxknots.v)>1){
    stop("Error: maxknots.m or maxknots.v greater than length 1")
  }
  control<-do.call(VarReg.control, control)
  n<-length(y)
  m<-list()
  col.n<-paste("Mean", c("zero", "constant", "linear", paste("Knot",seq(0,maxknots.m), sep="")), sep="_")
  row.n<-paste("Var", c("constant", "linear", paste("Knot",seq(0,maxknots.v), sep="")), sep="_")
  aic<-matrix(NA, ncol=maxknots.m+4,nrow=maxknots.v+3, dimnames=list(row.n, col.n))
  aicc<-matrix(NA, ncol=maxknots.m+4,nrow=maxknots.v+3, dimnames=list(row.n, col.n))
  bic<-matrix(NA, ncol=maxknots.m+4,nrow=maxknots.v+3, dimnames=list(row.n, col.n))
  hqc<-matrix(NA, ncol=maxknots.m+4,nrow=maxknots.v+3, dimnames=list(row.n, col.n))
  ll<-matrix(NA, ncol=maxknots.m+4,nrow=maxknots.v+3, dimnames=list(row.n, col.n))

  rep<-0
  for (i in 1:(maxknots.m+4)){
    if (i==1){
      meanmodel<-"zero"
      knots.m=NULL
    }else if (i==2){
      meanmodel<-"constant"
    }else if (i==3){
      meanmodel<-"linear"
    }else if (i>3){
      meanmodel<-"semi"
      knots.m<-i-4
    }
    for (j in 1:(maxknots.v+3)){
      if (j==1){
        varmodel<-"constant"
        knots.v=NULL
      }else if (j==2){
        varmodel<-"linear"
      }else if (j>2){
        varmodel<-"semi"
        knots.v<-j-3
      }
      rep<-rep+1
      #print(paste("rep", rep))
      #print(paste("mean =", i, "   var=", j ))
      param<-i-1+j
      #print(paste("parameters", param))
      log <- capture.output({
      m[[rep]]<-semiVarReg(y=y,x=x, meanmodel=meanmodel, varmodel=varmodel,knots.m=knots.m, knots.v=knots.v,degree=degree, mono.var = mono.var, control=control)
      })
      m[[rep]]$pos<-c(j,i)
      ll[j,i]<-m[[rep]]$loglik
      aicc[j,i]<-(2*(param)*n)/(n-(param)-1)-2*m[[rep]]$loglik
      aic[j,i]<-(2*(param))-2*m[[rep]]$loglik
      bic[j,i]<-log(n)*param-2*m[[rep]]$loglik
      hqc[j,i]<-log(log(n))*param-2*m[[rep]]$loglik
    }
  }

  if (selection=="BIC"){
    t<-which(bic == min(bic, na.rm=TRUE), arr.ind=TRUE)
  }else if (selection=="AICc"){
    t<-which(aicc == min(aicc, na.rm=TRUE), arr.ind=TRUE)
  }else if (selection=="AIC"){
    t<-which(aic == min(aic, na.rm=TRUE), arr.ind=TRUE)
  }else if (selection=="HQC"){
    t<-which(hqc == min(hqc, na.rm=TRUE), arr.ind=TRUE)
  }

  best.model<-list()
  for (r in 1:rep){
    if (m[[r]]$pos[1]==t[1] & m[[r]]$pos[2]==t[2]){
      best.model<-m[[r]]
      break
    }
  }

  list(ll=ll, AIC=aic, AICc=aicc,BIC=bic,HQC=hqc, best.model=best.model)
}