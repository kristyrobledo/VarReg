#' Censored Linear mean and variance regression
#'
#' \code{censlinVarReg} performs censored multivariate mean and multivariate variance regression.
#' This function is designed to be used by the  \code{\link{semiVarReg}} function.
#' @param dat Dataframe containing outcome and covariate data. Outcome data must be in the first column, with censored values set to the limits. Covariates for mean and variance model in next columns.
#' @param var.ind Vector containing the column numbers of the data in 'dat' to be fit as covariates in the variance model. FALSE indicates constant variance option.
#' @param mean.ind Vector containing the column numbers of the data in 'dat' to be fit as covariates in the mean model. 0 indicates constant mean option. NULL indicates zero mean option.
#' @param cens.ind Vector containing the column number of the data in 'dat' to indicate the censored data. 0 indicates no censoring, -1 indicates left (lower) censoring and 1 indicates right (upper) censoring.
#' @param mean.intercept Logical to indicate if an intercept is to be included in the mean model. Default is TRUE.
#' @param para.space Parameter space to search for variance parameter estimates. "positive" means only search positive parameter space, "negative" means search only negative parameter space and "all" means search all. Default is all.
#' @param mean.init Vector of initial estimates to be used for the mean model.
#' @param var.init Vector of initial estimates to be used for the variance model.
#' @param control List of control parameters. See \code{\link{VarReg.control}}.
#' @param ... arguments to be used to form the default control argument if it is not supplied
#' directly
#' @return \code{censlinVarReg} returns a list of output including:
#' \itemize{
#' \item\code{converged}: Logical argument indicating if convergence occurred.
#' \item\code{iterations}: Total iterations performed of the EM algorithm.
#'  \item\code{reldiff}: the positive convergence tolerance that occured at the final iteration.
#'  \item\code{loglik}: Numeric variable of the maximised log-likelihood.
#'  \item\code{boundary}: Logical argument indicating if estimates are on the boundary.
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
#'  \item\code{data}: Dataframe containing the variables included in the model.
#'  }
#'
#'@export


censlinVarReg<-function(dat, mean.ind=c(2), var.ind=c(2), cens.ind=c(3), mean.intercept=TRUE, para.space=c("all", "positive", "negative"),mean.init=NULL, var.init=NULL, control=list(...), ...){
  para.space<-match.arg(para.space)
  control<-do.call(VarReg.control, control)
  loops<-list()
  ll<-vector()
  X<-dat[,1]
  censor.ind<-dat[,cens.ind]
  n<-length(X)
  totalx<-var.ind
  if (is.null(mean.ind[1])==TRUE){
    meanmodel<-NULL
  }else if(mean.ind[1]==0){
    meanmodel<-FALSE
  }else if (mean.ind[1]>0){
    meanmodel<-data.frame(dat[,mean.ind])
  }

  ### constant variance first
  if (var.ind[1]==FALSE){
    loops[[1]]<-censloop_em(meanmodel, theta.old=1,beta.old=1, p.old=rep(1, n), x.0=NULL, censor.ind=censor.ind, X=X,mean.intercept=mean.intercept, maxit=control$maxit, eps=control$eps)
     var<- rep(loops[[1]]$theta.new, n)
    ll<- -n/2*log(2*pi)-1/2*sum(log(var))-sum(((X-loops[[1]]$fitted)**2)/(2*(var)))
    if (loops[[1]]$conv==FALSE){
      writeLines(paste("Warning: Did not converge at Maxit=", control$maxit))
    }
    name_alpha<-c("Intercept")
  ##variance model option
  }else{
    ## set parameters
    minmax<-list()
    x<-as.matrix(dat[,var.ind])
    if (is.null(var.init)==TRUE) {
      theta.old<-rep(1, 1+length(var.ind))
    }else if (is.null(var.init)==FALSE) {
      if (length(var.init)!=(1+length(var.ind))){
        stop("Error: length of var.init doesnt match model requirements. (length of vector not correct for # of parameters)")
      } else {
        theta.old<-var.init
      }
    }
    if (is.null(mean.init)==TRUE){
      beta.old<-rep(1, 1+length(mean.ind))
    }else if (is.null(mean.init)==FALSE){
      if (mean.intercept==TRUE && length(mean.init)!=(1+length(mean.ind))) {
        stop("Error: length of mean.init doesnt match model requirements. (length of vector not correct for # of parameters)")
      } else if (mean.intercept==FALSE && length(mean.init)!=(length(mean.ind))){
        stop("Error: length of mean.init doesnt match model requirements. (length of vector not correct for # of parameters). Note you are fitting a mean.intercept=FALSE model")
        } else {
        beta.old<-mean.init
      }
    }

    minmax[[1]]<-c("Min", "Max")
    ##find all combinations that need to be performed

    if (para.space=="all"){
      comb<-expand.grid(rep(minmax, times=length(var.ind)))
    }else if(para.space=="negative"){
      comb<-as.data.frame(matrix(c(rep("Max", length(rep))), 1, length(totalx), byrow=TRUE))
    }else if(para.space=="positive"){
      comb<-as.data.frame(matrix(c(rep("Min", length(rep))), 1, length(totalx), byrow=TRUE))
    }


    for (i in 1:nrow(comb)){
      p.old<-rep(theta.old[1], n)
      x.0<-matrix(NA, ncol=ncol(x), nrow=nrow(x) )
      for (j in 1:ncol(comb)){
        if(comb[i,j]=="Min"){
          x.0[,j]<-x[,j]-min(x[,j])
          p.old<-p.old+(theta.old[j]*x.0[,j])
        }
        if(comb[i,j]=="Max"){
          x.0[,j]<-max(x[,j])-x[,j]
          p.old<-p.old+(theta.old[j]*x.0[,j])
        }
      }
      loops[[i]]<-censloop_em(meanmodel, theta.old,beta.old, p.old, x.0, X,censor.ind, mean.intercept, control$maxit, control$eps)

      ll[i]<- -n/2*log(2*pi)-1/2*sum(log(loops[[i]]$p.old))-sum(((X-loops[[i]]$fittedmean)**2)/(2*(loops[[i]]$p.old)))

      n0=length(censor.ind[censor.ind==0])
      ll_1<-vector()
      ll1<-vector()
      p1<-vector()
      p2<-vector()
      p1[censor.ind==0]<- (-1/2*(log(loops[[i]]$p.old)))[censor.ind==0]
      p2[censor.ind==0]<- (-1*(((X-loops[[i]]$fittedmean)**2)/(2*(loops[[i]]$p.old))) )[censor.ind==0]
      ll0<- (-n0/2*log(2*pi)-1/2*sum(log(loops[[i]]$p.old)) -sum(((X-loops[[i]]$fittedmean)**2)/(2*(p.old))) )[censor.ind==0]
      ll0<- (-n0/2*log(2*pi)+sum(p1,p2,na.rm=TRUE))

      ll_1[censor.ind==-1]<- (log(stats::pnorm((X-loops[[i]]$fittedmean)/sqrt(loops[[i]]$p.old))))[censor.ind==-1]
      ll1[censor.ind==1]<- (log(1-stats::pnorm((X-loops[[i]]$fittedmean)/sqrt(loops[[i]]$p.old))))[censor.ind==1]
      ll[i]<-sum(ll0,ll_1,ll1, na.rm=TRUE)

      if (loops[[i]]$conv==FALSE){
        writeLines(paste("Warning: Did not converge at maxit=", control$maxit))
      }
    }
    name_alpha<-c("Intercept", colnames(dat)[var.ind])
  }
    ##find highest LL for both constant and nonconstant!
    max.ll<-which.max(ll)
    alpha<-vector()
    ##save all estimates
    alpha[1]<- loops[[max.ll]]$theta.new[1]
    conv.final<-loops[[max.ll]]$conv
    it.final<-loops[[max.ll]]$it
    reldiff.final<-loops[[max.ll]]$reldiff
    ll.final<-ll[[max.ll]]
    mean<-loops[[max.ll]]$mean
    if (is.na(reldiff.final)==TRUE){
      fit<- list(converged=conv.final,iterations=it.final,reldiff=reldiff.final, loglik=ll.final, boundary=NULL, aic.c=NULL, aic=NULL,bic=NULL,hqc=NULL,mean.ind=mean.ind,mean=NULL,var.ind=var.ind, variance=NULL,cens.ind=cens.ind, data=dat)
      return(fit)
    }

   if (var.ind[1]!=FALSE){
     for (j in 1:ncol(comb)){
      if(comb[max.ll,j]=="Min"){
        alpha[1]<-alpha[1]-loops[[max.ll]]$theta.new[j+1]*min(x[,j])
        alpha[j+1]<- loops[[max.ll]]$theta.new[j+1]
      }
      if(comb[max.ll,j]=="Max"){
        alpha[1]<-alpha[1]+loops[[max.ll]]$theta.new[j+1]*max(x[,j])
        alpha[j+1]<- -loops[[max.ll]]$theta.new[j+1]
      }
     }
   }
    names(alpha)<-name_alpha

   if (is.null(mean.ind)==FALSE){
     if (mean.intercept==TRUE){
       names(mean)<-c("Intercept", colnames(dat)[mean.ind])
     }else if (mean.intercept==FALSE){
       names(mean)<-c(colnames(dat)[mean.ind])
     }
   }
  if (is.null(mean.ind[1])==TRUE){
    param<-length(totalx)+1
  }else if (mean.ind[1]==0){
    param<-1+length(totalx)+1
  }else{
    if (mean.intercept==TRUE){
      param<-length(mean.ind)+1+length(totalx)+1
    }else if (mean.intercept==FALSE){
      param<-length(mean.ind)+length(totalx)+1
    }
  }

  if (sum(as.integer(loops[[max.ll]]$p.old < control$bound.tol))>0){
    boundary=TRUE
  }else{
    boundary=FALSE
  }
    ic<-criterion(n, ll.final, param)

  fit<- list(converged=conv.final,iterations=it.final,reldiff=reldiff.final, loglik=ll.final, boundary=boundary, aic.c=ic$aicc, aic=ic$aic,bic=ic$bic,hqc=ic$hqc,mean.ind=mean.ind,mean=mean,var.ind=var.ind, variance=alpha,cens.ind=cens.ind, data=dat)
  return(fit)
}
