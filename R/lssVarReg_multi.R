#' Semi parametric location, shape and scale regression
#'
#' \code{lssVarReg.multi} performs a semiparametric location (\eqn{\xi} or xi), shape (\eqn{\nu} or nu) and scale (\eqn{\omega} or omega) regression model. This is designed for multiple covariates that are fit in the location, scale and shape models.
#' @param y Vector containing outcome data. Must be no missing data.
#' @param x Vector containing the covariate data, same length as \code{y}. Must be no missing data.
#' @param locationmodel Text to specify the location model to be fit. Options: \code{"constant"} = constant model (intercept only), \code{"linear"} = linear term with x covariate, \code{"semi"} = semiparametric spline (specify with \code{knots.l}).
#' @param scale2model Text to specify the scale^2 model to be fit. Options: \code{"constant"} = constant term only, \code{"linear"} = linear term with \code{x} covariate, \code{"semi"} = semiparametric spline (specify with \code{knots.sc})
#' @param shapemodel Text to specify the shape model to be fit. Options: \code{"constant"} = constant shape model, \code{"linear"} = linear term with x covariate, \code{"semi"} = semiparametric spline (specify with \code{knots.sh}).
#' @param knots.l Integer indicating the number of internal knots to be fit in the location model. Default is '2'. (Note that the knots are placed equidistantly over x.)
#' @param knots.sc Integer indicating the number of internal knots to be fit in the scale^2 model. Default is '2'. (Note that the knots are placed equidistantly over x.)
#' @param knots.sh Integer indicating the number of internal knots to be fit in the shape model. Default is '2'. (Note that the knots are placed equidistantly over x.)
#' @param degree Integer to indicate the degree of the splines fit in the location and scale.  Default is '2'.
#' @param mono.scale Text to indicate whether the scale2 model is monotonic. Default is \code{"none"} (no monotonic constraints). Options are \code{"inc"} for increasing or \code{"dec"} for decreasing. If this is chosen, the appropriate \code{para.space} is set automatically (\code{"positive"} for \code{inc}, \code{"negative"} for \code{dec}).
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


lssVarReg.multi<-function(y, x,
                          locationmodel=c("constant", "linear", "semi"),
                          location.vars = c(1),
                          scale2model=c("constant", "linear", "semi"),
                          scale2.vars = c(1),
                          shapemodel=c("constant", "linear", "semi"),
                          shape.vars = c(1),
                          knots.l=NULL, knots.sc=NULL, knots.sh=NULL, degree=2,
                          location.init=NULL, scale2.init=NULL,shape.init=NULL,
                          int.maxit=1000, print.it=FALSE, control=list(...), ...) {
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
  if (locationmodel[1]=="constant"){
    knots.l<-NULL
    location.vars<-NULL
  } else if (length(locationmodel)!=length(location.vars) ||
      sum(locationmodel=="semi")!=length(knots.l) ){
    stop("Error: locationmodel, knots.l and location.vars not equal length")
  }
  if (sum(location.vars=="semi")>0 &&
      sum(location.vars=="semi")!=length(knots.l)){
    stop("Error: vector knots.l not long enough")
  }

  if (scale2model[1]=="constant"){
    knots.sc=NULL #ensure these are set correctly
    scale2.vars=NULL
  } else if (length(scale2model)!=length(scale2.vars) ||
             sum(scale2model=="semi")!=length(knots.sc)  ){
    stop("Error: scale2model, knots.sc and scale2.vars not equal length")
  }
  if (sum(scale2model=="semi")>0 &&
      sum(scale2model=="semi")!=length(knots.sc))(
        stop("Error: vector knots.sc not long enough")
      )

  if (shapemodel[1]=="constant"){
    knots.sh=NULL #ensure these are set correctly
    shape.vars=NULL
  } else if (length(shapemodel)!=length(shape.vars) ||
             sum(shapemodel=="semi")!=length(knots.sh)  ){
    stop("Error: shapemodel, knots.sh and shape.vars not equal length")
  }
  if (sum(shapemodel=="semi")>0 &&
      sum(shapemodel=="semi")!=length(knots.sh))(
        stop("Error: vector knots.sh not long enough")
      )


  if ("semi" %in% locationmodel && length(which("semi"== locationmodel))!=length(knots.l)){
    stop("Error: please specify the number of knots for each semiparametric variable in the location model")
  }
  if ("semi" %in% scale2model && length(which("semi"== scale2model))!=length(knots.sc)){
    stop("Error: please specify the number of knots for each semiparametric variable in the scale2 model")
  }
  if ("semi" %in% shapemodel && length(which("semi"== shapemodel))!=length(knots.sh)){
    stop("Error: please specify the number of knots for each semiparametric variable in the shape model")
  }
  n<-length(y)

  alldat<-data.frame(y, mean.int=rep(1,n), x)
  colnames(x)<-make.names(colnames(x))

    #loop thru the location variables
  mean.ind<-NULL
  msemicounter<-0
  if (locationmodel[1]=="constant"){
    mean.ind<-c(2)
    xiold<-1
  }else if (length(locationmodel)>=1){
    mean.ind<-c(2)
    for (i in 1:length(locationmodel)){
      #print(i)
      if(locationmodel[i]=="semi"){
        print("semi")
        msemicounter<-msemicounter+1
        bmean<-splines::bs(x=x[,location.vars[i]], df=(degree+knots.l[msemicounter]), degree=degree)
        colnames(bmean) <- paste(paste(paste(colnames(x)[location.vars[i]], "Knt", sep="_"),knots.l[msemicounter], sep = ""), paste("Base", colnames(bmean), sep=""), sep = "_")
        alldat<-data.frame(alldat, bmean)
        mean.ind[length(mean.ind)+1:(ncol(bmean))]<-which(colnames(alldat)%in%colnames(bmean))
       }else if (locationmodel[i]=="linear"){
        print("linear")
        mean.ind[length(mean.ind)+1]<-location.vars[i]+2 ##assign next free place in vector with covariate
        }
    }
    if (is.null(location.init)==TRUE){
      xiold<-rep(1,times = length(mean.ind))
    }else if (length(location.init)==length(mean.ind)){
      xiold<-location.init
    }else{
            stop("Error: check location.init is the correct length (expecting intercept + parameter starting estimates)")
     }
     }else {stop("Error: check locationmodel contains the appropriate strings")
    }

  #loop thru the scale2 variables
  mono.scale<-"none"
  var.ind<-NULL
  vsemicounter<-0
  if (scale2model[1]=="constant"){
    knots.sc<-NULL
    var.ind<-FALSE
    omega2old<-1
    }else if (length(scale2model)>=1){
      var.ind<-NULL

      for (i in 1:length(scale2model)){
        print(i)
        if(scale2model[i]=="semi"){
          print("semi")
          vsemicounter<-vsemicounter+1
          bvar<-splines::bs(x=x[,scale2.vars[i]], df=(degree+knots.sc[vsemicounter]), degree=degree)
          colnames(bvar) <- paste(paste(paste(colnames(x)[scale2.vars[i]],"Knt",sep="_"),knots.sc[vsemicounter], sep = ""), paste("Base", colnames(bvar), sep=""), sep = "_")
          alldat<-data.frame(alldat, bvar)
          var.ind[length(var.ind)+1:(ncol(bvar))]<-which(colnames(alldat)%in%colnames(bvar))
         }else if (scale2model[i]=="linear"){
          print("linear")
          var.ind[length(var.ind)+1]<-scale2.vars[i]+2 ##assign next free place in vector with covariate
         }
      }
      if (is.null(scale2.init)==TRUE){
        omega2old<-rep(1,times = 1+length(var.ind))
      }else if (length(scale2.init)==1+length(var.ind)){
        omega2old<-scale2.init
      }else{
        stop("Error: check scale2.init is the correct length (expecting intercept + parameter starting estimates)")
      }


    }else {stop("Error: check scale2model contains the appropriate strings")
    }


    ssemicounter<-0
    nu.ind<-NULL
    if (shapemodel[1]=="constant"){
    nuold<-1
    knots.sh<-NULL
    nu.ind<-NULL
  } else if (length(shapemodel)>=1){
    nu.ind<-NULL
    for (i in 1:length(shapemodel)){
      print(i)
      if(shapemodel[i]=="semi"){
        print("semi")
        ssemicounter<-ssemicounter+1
        bsh<-splines::bs(x=x[,shape.vars[i]], df=(degree+knots.sh[ssemicounter]), degree=degree)
        colnames(bsh) <- paste(paste(paste(colnames(x)[shape.vars[i]],"Knt",sep="_"),knots.sh[ssemicounter], sep = ""), paste("Base", colnames(bsh), sep=""), sep = "_")
        alldat<-data.frame(alldat, bsh)
        nu.ind[length(nu.ind)+1:(ncol(bsh))]<-which(colnames(alldat)%in%colnames(bsh))
      }else if (shapemodel[i]=="linear"){
        nu.ind[length(nu.ind)+1]<-shape.vars[i]+2 ##assign next free place in vector with covariate
      }
    }
    nuold<-rep(1,times = 1+length(nu.ind))
    if (is.null(shape.init)==TRUE){
      nuold<-rep(1,times = 1+length(nu.ind))
    }else if (length(shape.init)==1+length(nu.ind)){
      nuold<-shape.init
    }else{
      stop("Error: check shape.init is the correct length (expecting intercept + parameter starting estimates)")
    }
  }else {stop("Error: check shapemodel contains the appropriate strings")
  }


    l<-loop_lss(alldat,xiold,omega2old,nuold,mean.ind, var.ind, nu.ind, para.space="all",
                maxit=control$maxit, eps=control$epsilon,int.maxit, print.it)
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
    if (shapemodel[1]=="constant"){
      names(l$nunew)<-"Intercept"
    } else {
      names(l$nunew)<-c("Intercept", colnames(alldat)[nu.ind])
    }

  out<-list(modeltype="LSS model", locationmodel=locationmodel, knots.l=knots.l, scale2model=scale2model, knots.sc=knots.sc, shapemodel=shapemodel, knots.sh=knots.sh, degree=degree, converged=l$conv, iterations=l$it,reldiff=l$reldiff, loglik=loglik, aic.c=ic$aicc, aic=ic$aic,bic=ic$bic, mono.scale=mono.scale, hqc=ic$hqc,
            location=l$xinew,
            scale2=l$omega2new,
            shape=l$nunew, data=alldat)
  class(out) <- c("lssVarReg")
  return(out)
}

