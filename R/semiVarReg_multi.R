#' Semi parametric mean and variance regression (multivariate)
#'
#' \code{semiVarReg.multi} performs semi-parametric mean and variance regression models. This is
#' designed for multiple covariates fit in the mean and variance models.
#'@param y Vector containing outcome data. Must be no missing data and any censored values must
#' be set to the limits of detection.
#'@param x Matrix containing the covariate data. Must be no missing data and same length as \code{y}.
#'@param meanmodel Vector to specify the mean model to be fit to the data. The possible inputs are
#'  \code{"zero"}, \code{"constant"}, or a vector to indicate if covariates are to be \code{"linear"} or \code{"semi"}. \code{"semi"}
#'  indicates a semi-parametric spline model, with the number of internal knots specified in
#'  \code{knots.m}.
#'@param varmodel Vector to specify the variance model to be fit to the data. The possible inputs are
#' \code{"constant"}, or a vector to indicate if each covariate is to be \code{"linear"} or \code{"semi"}. \code{"semi"} indicates a semi-parametric
#' B-spline model, with the number of internal knots specified in \code{knots.v}.
#'@param knots.m Integer indicating the number of internal knots to be fit in the semi-parametric
#' mean model. Knots are placed equidistantly over the covariate. The default value is \code{2}.
#'@param knots.v Integer indicating the number of internal knots to be fit in the semi-parametric
#'  variance model. Knots are placed equidistantly over the covariate. The default value is \code{2}.
#'@param degree Integer indicating the degree of the splines fit in the mean and the variance models.
#' The default value is \code{2}.
#' @param control list of control parameters. See \code{\link{VarReg.control}}.
#' @param ... arguments to be used to form the default control argument if it is not supplied
#' directly
#' @return \code{semiVarReg.multi} returns an object of class \code{"VarReg"} which inherits some components from the class \code{"glm"}. This object of class \code{"VarReg"} is a list containing the following components:
#' \itemize{
#'  \item\code{modeltype}: Text indicating the model that was fit, indicating an uncensored approach was performed.
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
#'  \item\code{data}: Dataframe containing the variables included in the model.}
#'
#'@examples
#'data(mcycle)
#'## run a model with linear mean and linear variance:
#'linmodel<-semiVarReg.multi(mcycle$accel, mcycle$times, meanmodel="linear", varmodel="linear",
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


semiVarReg.multi<-function(y,x,
                           mean.model=c("zero", "constant", "linear", "semi"),
                           mean.vars=c(1),knots.m=NULL,
                           var.model=c("constant", "linear", "semi"),
                           var.vars=c(1),
                           knots.v=NULL, degree=2,
                           control=list(...), ...){
  control<-do.call(VarReg.control, control)
  if (anyNA(x)){
    stop("Error: x contains missing data")
  }
  if (anyNA(y)){
    stop("Error: y contains missing data")
  }
  if (length(y)!=nrow(x)){
    stop("Error: x and y not same length")
  }
  #check lengths of all mean model components are equal
  if ((mean.model=="zero" || mean.model=="constant")){
    knots.m=NULL
    mean.vars=NULL
  } else if (length(mean.model)!=length(mean.vars) ||
             sum(mean.model=="semi")!=length(knots.m)){
    stop("Error: mean.model, knots.m and mean.vars not equal length")
  }
  if (sum(mean.model=="semi")>0 &&
      sum(mean.model=="semi")!=length(knots.m))(
        stop("Error: vector knots.m not long enough")
      )

  if (var.model=="zero" || var.model=="constant"){
    knots.v=NULL #ensure these are set correctly
    var.vars=NULL
  } else if (length(var.model)!=length(var.vars) ||
      sum(var.model=="semi")!=length(knots.v)  ){
    stop("Error: var.model, knots.v and var.vars not equal length")
  }
  if (sum(var.model=="semi")>0 &&
      sum(var.model=="semi")!=length(knots.v))(
        stop("Error: vector knots.v not long enough")
      )


  if ("semi" %in% mean.model && length(which("semi"== mean.model))!=length(knots.m)){
    stop("Error: please specify the number of knots for each semiparametric variable in the mean model")
    }
  if ("semi" %in%var.model && length(which("semi"== var.model))!=length(knots.v)){
    stop("Error: please specify the number of knots for each semiparametric variable in the variance model")
    }

  alldat<-data.frame(y,x)
  colnames(x)<-make.names(colnames(x))

  #loop thru the mean variables
  mean.ind<-NULL
  msemicounter<-0
  if (mean.model[1]=="constant"){
    mean.ind<-0
    knots.m<-NULL
  }else if (mean.model[1]=="zero"){
    mean.ind<-NULL
    knots.m<-NULL
  } else if (length(mean.model)>=1){
    mean.ind<-NULL
      for (i in 1:length(mean.model)){
        #print(i)

        if(mean.model[i]=="semi"){
          #print("semi")
          msemicounter<-msemicounter+1
          bmean<-splines::bs(x=x[,mean.vars[i]], df=(degree+knots.m[msemicounter]), degree=degree)
          colnames(bmean) <- paste(paste(paste(colnames(x)[mean.vars[i]], "Knt",sep="_"),knots.m[msemicounter], sep = ""), paste("Base", colnames(bmean), sep=""), sep = "_")
          alldat<-data.frame(alldat, bmean)
          mean.ind[length(mean.ind)+1:(ncol(bmean))]<-which(colnames(alldat)%in%colnames(bmean))
        }else if (mean.model[i]=="linear"){
          #print("linear")
          mean.ind[length(mean.ind)+1]<-mean.vars[i]+1 ##assign next free place in vector with covariates
        }
      }
  }else {stop("Error: check mean.model contains the appropriate strings")
  }



  #loop thru the variance variables
  var.ind<-NULL
  vsemicounter<-0
  if (var.model[1]=="constant"){
    var.ind<-FALSE
    knots.v<-NULL
   } else if (length(var.model)>=1){
    var.ind<-NULL
    for (i in 1:length(var.model)){
     # print(i)
      if(var.model[i]=="semi"){
      #  print("semi")
        vsemicounter<-vsemicounter+1
        bvar<-splines::bs(x=x[,var.vars[i]], df=(degree+knots.v[vsemicounter]), degree=degree)
        colnames(bvar) <- paste(paste(paste(colnames(x)[var.vars[i]], "Knt", sep="_"),knots.v[vsemicounter], sep = ""), paste("Base", colnames(bvar), sep=""), sep = "_")
        alldat<-data.frame(alldat, bvar)
        var.ind[length(var.ind)+1:(ncol(bvar))]<-which(colnames(alldat)%in%colnames(bvar))
      }else if (var.model[i]=="linear"){
       # print("linear")
        var.ind[length(var.ind)+1]<-var.vars[i]+1 ##assign next free place in vector with covariate
      }
    }
  }else {stop("Error: check var.model contains the appropriate strings")
  }


#run the model by calling linVarReg function
 # if (is.null(cens.ind)==TRUE){
    result<-linVarReg(alldat, var.ind=var.ind, mean.ind=mean.ind, para.space="all", control=control)
    type<-"Mean and Variance regression"
#  }else {
#    alldat<-cbind(alldat,cens.ind)
#    c.ind<-which(colnames(alldat)%in%("cens.ind"))
#    result<-censlinVarReg(dat=alldat, var.ind=var.ind, mean.ind=mean.ind,mean.intercept=mean.intercept, cens.ind=c.ind, para.space=para.space, control=control)
#    type<-c("Censored Mean and Variance regression with mean.intercept=",mean.intercept)
#  }

  model<-list(modeltype=type, knots.m=knots.m, knots.v=knots.v, degree=degree,
              meanmodel=mean.model, varmodel=var.model)
  result<-c(model, result)
  class(result) <- c("VarReg")
  return(result)
}
