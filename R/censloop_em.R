#' The Censored data EM loop
#'
#' \code{censloop_em} is an EM loop function for censored data to be utilised by various other higher level functions.
#' @param meanmodel Dataframe containing only the covariates to be fit in the mean model. NULL for zero mean model and FALSE for constant mean model.
#' @param theta.old Vector containing the initial variance parameter estimates to be fit in the variance model.
#' @param beta.old Vector containing the initial mean parameter estimates to be fit in the mean model.
#' @param p.old Vector of length n containing the initial variance estimate.
#' @param x.0 Matrix of covariates (length n) to be fit in the variance model. All have been rescaled so zero is the minimum. If NULL, then its a constant variance model.
#' @param X Vector of length n of the outcome variable.
#' @param censor.ind Vector of length n of the censoring indicator. 0=uncensored, -1=left censored and 1 is right censored.
#' @param mean.intercept Logical to indicate if mean intercept is to be included in the model.
#' @param maxit Number of maximum iterations for the EM algorithm.
#' @param eps Very small number for the convergence criteria.
#' @return A list of the results from the EM algorithm, including:
#' \itemize{
#' \item\code{conv}: Logical argument indicating if convergence occurred
#' \item\code{it}: Total iterations performed of the EM algorithm
#'  \item\code{reldiff}: the positive convergence tolerance that occured at the final iteration.
#'  \item\code{theta.new}: Vector of variance parameter estimates. Note that these are not yet
#'  transformed back to the appropriate scale
#'  \item\code{mean}: Vector of mean parameter estimates
#'  \item\code{fittedmean}: Vector of fitted mean estimates
#'  \item\code{p.old}: Vector of fitted variance estimates
#'  }
#' @export

censloop_em<-function(meanmodel, theta.old,beta.old, p.old, x.0, X,censor.ind, mean.intercept, maxit, eps){
  n<-length(X)
  R2<-list()
  Z20<-list()
  Z21<-list()
  Y2<-rep(NA, times=n)
  Ca<- rep(NULL, times=n)
  conv<-FALSE
  it<-NULL
  p.old[p.old < 1e-10] <- 1e-10 ##finite bound on the variance
  p.old[p.old > 1e10] <- 1e10 ##finite bound on the variance
  theta.new<-theta.old
  z=sqrt(1/p.old)
  c_1<-X*z
  c_1[censor.ind==-1] <- NA
  c_2<-X*z
  c_2[censor.ind==1]<-NA
  beta.old<-NULL
  outmean<-m_models(c_1, c_2, z,n, meanmodel,mean.intercept, beta.old, maxit, eps)
  mean.fit<-outmean$mean.fit
  beta.old<-outmean$beta.new
  beta.new<-beta.old
  for (q in 1:maxit) {
    it<-q
    X0<-X-mean.fit
    Y2[censor.ind==0]<-(theta.old[1]+(theta.old[1]**2/p.old)*(X0**2/p.old - 1))[censor.ind==0]
    Ca[censor.ind==-1 & round(stats::pnorm((X0)/sqrt(p.old)), digits=5)!=0]<-((stats::pnorm((X0)/sqrt(p.old)))/(-stats::dnorm((X0)/sqrt(p.old))))[censor.ind==-1 & round(stats::pnorm((X0)/sqrt(p.old)), digits=5)!=0]
    Ca[censor.ind==-1 & round(stats::pnorm((X0)/sqrt(p.old)), digits=5)==0]<- 1/((X0)/sqrt(p.old))[censor.ind==-1 & round(stats::pnorm((X0)/sqrt(p.old)), digits=5)==0] ##accurate numerical representation when x tends to zero
    Y2[censor.ind==-1]<-(theta.old[1]+(theta.old[1]**2/p.old)*((X0/sqrt(p.old))/Ca))[censor.ind==-1]
    Ca[censor.ind==1 & round(stats::pnorm((X0)/sqrt(p.old)), digits=5)!=0]<-((1-stats::pnorm((X0)/sqrt(p.old)))/(stats::dnorm((X0)/sqrt(p.old))))[censor.ind==1 & round(stats::pnorm((X0)/sqrt(p.old)), digits=5)!=0]
    Ca[censor.ind==1 & round(stats::pnorm((X0)/sqrt(p.old)), digits=5)==0]<- 1/((X0)/sqrt(p.old))[censor.ind==1 & round(stats::pnorm((X0)/sqrt(p.old)), digits=5)==0]
    Y2[censor.ind==1]<- (theta.old[1]+(theta.old[1]**2/p.old)*((X0/sqrt(p.old))/Ca))[censor.ind==1]

  theta.new[1]<-mean(Y2)

  if (!is.null(x.0)){
    X00<-X0[censor.ind==0]
    X01<-X0[censor.ind!=0]
    x00<-x.0[censor.ind==0,]
    x01<-x.0[censor.ind!=0,]
    if (length(theta.old)<3){
      x00<-as.matrix(x00)
      x01<-as.matrix(x01)
    }
    pold0<-p.old[censor.ind==0]
    pold1<-p.old[censor.ind!=0]
    Ca1<-Ca[censor.ind!=0]
      Z20<-sapply(1:ncol(x00),function(x) x00[,x]*theta.old[x+1]+((theta.old[x+1]*x00[,x])**2/pold0)*(X00**2/pold0-1))
      Z21<-sapply(1:ncol(x01),function(x) x01[,x]*theta.old[x+1]+((theta.old[x+1]*x01[,x])**2/pold1)*(X01/sqrt(pold1))/Ca1)
      zz<-rbind(Z20, Z21)
      xx<-rbind(x00, x01)
      theta.new[-1]<-sapply(1:ncol(x.0), function(x) mean(zz[,x]/xx[,x],na.rm=TRUE ))
    }

    ##calculate weighting for next LM - bounded at 1e10
   z=1/sqrt(p.old)
   z[z > 1e10] <- 1e10 ##finite bound
   c_1<-X*z
   c_1[censor.ind==-1]<-NA
   c_2<-X*z
   c_2[censor.ind==1]<-NA

   outmean<-tryCatch({
     m_models(c_1, c_2, z,n, meanmodel,mean.intercept,beta.old, maxit, eps)
   },error=function(cond){
     message("Mean model (survreg) failed to converge at maxit=", maxit, ". Attempting more iterations.")
     return(NULL)
   },warning=function(cond){
     message("Mean model (survreg) failed to converge at maxit=", maxit, ". Attempting more iterations.")
   })
   if (is.null(outmean)==TRUE || any(is.nan(outmean$beta.new))){
     outmean<-tryCatch({
       m_models(c_1, c_2, z,n, meanmodel,mean.intercept,beta.old, maxit*100, eps)
     },error=function(cond){
       message("Mean model (survreg) still failed to converge at maxit=", maxit, "*100. Review initial estimates.")
       return(NULL)
     },warning=function(cond){
       message("Mean model (survreg) still failed to converge at maxit=", maxit, "*100. Review initial estimates.")
     })
   }
   if (is.null(outmean)==TRUE){
     break}
   beta.new<-outmean$beta.new
   mean.fit<-outmean$mean.fit
   if (all(is.finite(c(beta.new, theta.new)))==TRUE){
     reldiff<-sqrt(sum((c(theta.new,beta.new)-c(theta.old,beta.old))**2)/sum(c(beta.old,theta.old))**2)
   }else {
     reldiff<-NaN
     print("Estimates are not finite. Please review initial estimates.")
     print(c("Mean:",beta.new))
     print(c("Variance",theta.new))
     break
   }


   theta.old<-theta.new
   beta.old<-beta.new
   p.old<-rep(theta.old[1], n)
    if (!is.null(x.0)){
      p.old<-rowSums(cbind(rep(theta.new[1], n),sapply(1:ncol(x.0), function(x) theta.new[x+1]*x.0[,x])))
    }
   if (sum(p.old<0)>0){
     print("Estimates for variance/scale have gone negative. Please review initial estimates.")
    break
   }
   #print(p.old)
   #print(it)
   # print(c(beta.new, theta.new))
    if (is.finite(reldiff)==TRUE){
      if (reldiff<eps){
        conv<-TRUE
        break
      }
    }
  }
  if (is.null(outmean)==TRUE){
    list(conv=FALSE, reldiff=NA, it=it, mean=beta.old, theta.new=theta.old, fittedmean=NULL, p.old=p.old)
  }else{
  list(conv=conv, reldiff=reldiff, it=it, mean=beta.new, theta.new=theta.new, fittedmean=mean.fit, p.old=p.old)
  }
}


m_models<-function(c_1, c_2, z,n, meanmodel,mean.intercept,beta.old, maxit, eps){
  if (is.null(meanmodel)==TRUE){
    mean.fit<-rep(0,n)
    beta.new<-NULL
    }else if (is.data.frame(meanmodel)==TRUE){
      if (mean.intercept==TRUE){
      meanmodel2<-as.data.frame(apply(meanmodel, 2, function(x) x*z))
      l<-survival::survreg(survival::Surv(c_1, c_2, type="interval2")~ -1+z+ . , data=meanmodel2,  dist = "gaussian",  maxiter=maxit, rel.tolerance=eps, init=beta.old)
      beta.new<-l$coeff
      mean.fit<-l$linear.predictor/z
      }else if (mean.intercept==FALSE){
      meanmodel2<-as.data.frame(apply(meanmodel, 2, function(x) x*z))
      l<-survival::survreg(survival::Surv(c_1, c_2, type="interval2")~ -1+., data=meanmodel2,  dist = "gaussian", maxiter=maxit, rel.tolerance=eps, init=beta.old)
      beta.new<-l$coeff
      mean.fit<-l$linear.predictor/z
      }
      }else if (meanmodel[1]==FALSE){
    l<-survival::survreg(survival::Surv(c_1, c_2, type="interval2")~ -1+z, dist = "gaussian", maxiter=maxit, rel.tolerance=eps)
    beta.new<-l$coeff
    mean.fit<-l$linear.predictor/z
    }
  return(list(mean.fit=mean.fit, beta.new=beta.new))
  }