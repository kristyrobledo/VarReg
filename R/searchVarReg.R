#' Searches for best semi parametric mean and variance regression model
#'
#' \code{searchVarReg} performs multiple semi-parametric mean and variance regression models for a covariate of interest, in order to search for the optimal number of knots. The best model is chosen based on the information criterion of preference (\code{"selection"}). At the moment, this is only designed for a single covariate that is fit in both the mean and variance models.
#' @param y Vector containing outcome data. Must be no missing data and any censored values must
#' be set to the limits of detection.
#' @param x Vector containing the covariate data. Must be no missing data and same length as \code{y}.
#' @param cens.ind Vector containing the censoring indicator, if applicable. There must be no missing
#' data contained in the vector and this vector should be the same length as \code{y}.
#' \code{"0"} values indicate uncensored data, \code{"1"} indicates right, or upper, censoring and
#' \code{"-1"} indicates left, or lower, censoring. The default is \code{NULL} which indicates there
#' is no censored data.
#' @param maxknots.m  Integer indicating the maximum number of internal knots to be fit in the mean model. Default is \code{3}. (Note that the knots are placed equidistantly over x.)
#' @param maxknots.v Integer indicating the maximum number of internal knots to be fit in the variance model. Default is \code{3}. (Note that the knots are placed equidistantly over x.)
#' @param degree The degree of the splines fit in the mean and variance. Default is \code{2}.
#' @param mono.var Text to indicate whether the variance model is monotonic (only applied to 'linear' or
#'semi-parametric variance models). Default is "\code{none}" (no monotonic constraints). Options are
#'"\code{inc}" for increasing or "\code{dec}" for decreasing. If the variance model is linear, the
#'parameter space is constrained (positive for increasing and negative for decreasing). For
#'semi-parametric variance models, the appropriate monotonic B splines are fit in the
#'semi-parametric variance model.
#' @param selection Text to indicate which information criteria is to be used for the selection of the
#'  best model. Choices are "\code{AIC}", "\code{AICc}", "\code{BIC}" and "\code{HQC}".
#' Default is "\code{AIC}".
#' @param print.it Logical to indicate whether to print progress from each model as the models are
#'  performed. Default is \code{FALSE}.
#' @param control  list of control parameters. See \code{\link{VarReg.control}}.
#' @param ... arguments to be used to form the default control argument if it is not supplied
#' directly
#'
#' @return \code{searchVarReg} returns an list, with the following components:
#' \itemize{
#' \item \code{ll}: a dataframe of the log-likelihoods from each of the models that have been fit.
#' \item \code{AIC}: a dataframe of the AIC from each of the models that have been fit. The parameters
#'  fit in the mean model are given in the columns, and the parameters in the variance are given
#'  in the rows.
#'  \item \code{AICc}: a dataframe of the AIC-c from each of the models that have been fit.
#'  \item \code{BIC}: a dataframe of the BIC from each of the models that have been fit.
#'  \item \code{HQC}: a dataframe of the HQC from each of the models that have been fit.
#'  \item \code{best.model}: an object of class \code{VarReg} (see \code{\link{semiVarReg}})
#'  containing the output from the optimal model (that model within the specified models in
#'  the mean and variance with the lowest information criterion according to the criterion selected).
#' }
#'
#'
#'@details A matrix of models are performed, of increasing complexity. Mean models start at a zero mean
#' model, then constant mean, linear, 0 internal knots, etc, up to a maximum internal knots as specified
#'  in \code{maxknots.m}. Variance models start at constant variance, linear variance, 0 internal knots,
#'   etc, up to max internal knots as specified in \code{maxknots.v}.
#'
#' Note that this function can take some time to run, due to the number of models to be fit.
#' A window will appear on windows based systems to show a progress bar for the function.
#'
#'@seealso \code{\link{semiVarReg}}, \code{\link{VarReg.control}}
#'
#'@examples
#'data(mcycle)
#'### not run
#'### find<-searchVarReg(mcycle$accel, mcycle$times, maxknots.v=3, maxknots.m=3,
#'### selection="HQC", maxit=10000)

#'@export

searchVarReg<-function(y,x,cens.ind=NULL, maxknots.m=3, maxknots.v=3, degree=2, mono.var=c("none", "inc", "dec"), selection=c("AIC", "AICc", "HQC", "BIC"), print.it=FALSE, control=list(...), ...){
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
  pb <- utils::winProgressBar(title="Finding optimal model", label="0% done", min=0, max=100, initial=0)
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
      info <- sprintf("%d%% done", round((rep/((maxknots.m+4)*(maxknots.v+3)))*100))
      param<-i-1+j
      log <- utils::capture.output({
      m[[rep]]<-semiVarReg(y=y,x=x, cens.ind=cens.ind, meanmodel=meanmodel, varmodel=varmodel,knots.m=knots.m, knots.v=knots.v,degree=degree, mono.var = mono.var, control=control)
      })
      utils::setWinProgressBar(pb, rep/((maxknots.m+4)*(maxknots.v+3))*100, label=info)
      m[[rep]]$pos<-c(j,i)
      ll[j,i]<-m[[rep]]$loglik
      ic<-criterion(n=n, loglik=m[[rep]]$loglik, param=param)
      aicc[j,i]<-ic$aicc
      aic[j,i]<-ic$aic
      bic[j,i]<-ic$bic
      hqc[j,i]<-ic$hqc
      if (print.it==TRUE){
        print(c("Model", rep, " of ", sum((maxknots.m+4)*(maxknots.v+3))))
        print(c("Meanmodel= ",meanmodel))
        print(c("Varmodel=",  varmodel))
        print(c("LL=",  ll[j,i]))
        print(c("IC=",  ic))
        }
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
  close(pb)
  best.model<-list()
  for (r in 1:rep){
    if (m[[r]]$pos[1]==t[1] & m[[r]]$pos[2]==t[2]){
      best.model<-m[[r]]
      break
    }
  }

  list(ll=ll, AIC=aic, AICc=aicc,BIC=bic,HQC=hqc, best.model=best.model)
}