#' The EM loop
#'
#' \code{loop_em} is a basic EM loop function to be utilised by various other higher level functions
#' @param meanmodel Dataframe containing only the covariates to be fit in the mean model. NULL for zero mean model and FALSE for constant mean model.
#' @param theta.old Vector containing the initial variance parameter estimates to be fit in the variance model.
#' @param p.old Vector of length n containing the containing the initial variance estimate.
#' @param x.0 Matrix of covariates (length n) to be fit in the variance model. All have been rescaled so zero is the minimum. If NULL, then its a constant variance model.
#' @param X Vector of length n of the outcome variable.
#' @param maxit Number of maximum iterations for the EM algorithm, default 1000.
#' @param eps Very small number for the convergence criteria, default 1 times 10 power6
#' @return A list of the results from the EM algorithm, including conv, reldiff, it, mean, theta.new, fitted.
#'

loop_em<-function(meanmodel, theta.old, p.old, x.0, X, maxit, eps){
  n<-length(X)
  Z2<-list()
  R2<-list()
  conv<-FALSE
  it<-NULL
  theta.new<-theta.old
  if (is.null(meanmodel)==TRUE){
    mean.fit<-rep(0,n)
    beta.old<-NULL
  }else if (is.data.frame(meanmodel)==TRUE){
    l<-lm(X~., weight=1/p.old, data=meanmodel)
    beta.old<-l$coeff
    mean.fit<-l$fitted
  }else if (meanmodel[1]==FALSE){
    l<-lm(X~1, weight=1/p.old)
    beta.old<-l$coeff
    mean.fit<-l$fitted
  }
  beta.new<-beta.old
  for (q in 1:maxit) {
    it<-q
    X0<-X-mean.fit
    Y2<-theta.old[1]+(theta.old[1]**2/p.old)*(X0**2/p.old - 1)
    theta.new[1]<-mean(Y2)
    if (!is.null(x.0)){
      Z2<-sapply(1:ncol(x.0),function(x) x.0[,x]*theta.old[x+1]+((theta.old[x+1]*x.0[,x])**2/p.old)*(X0**2/p.old-1))
      theta.new[-1]<-sapply(1:ncol(x.0), function(x) mean(Z2[,x]/x.0[,x],na.rm=TRUE ))
    }

    ##calculate weighting for next LM - bounded at 1e10
    wt=1/p.old
    wt[wt > 1e10] <- 1e10
    ##calculate weighted mean - weighted by inverse variance FROM LM WEIGHTED MODEL
    if (is.null(meanmodel)==TRUE){
    }else if (is.data.frame(meanmodel)==TRUE){
      l<-lm(X~., weight=wt, data=meanmodel)
      beta.new<-l$coeff
      mean.fit<-l$fitted
    }else if (meanmodel[1]==FALSE){
      l<-lm(X~1, weight=wt)
      beta.new<-l$coeff
      mean.fit<-l$fitted
    }
    reldiff<-sqrt(sum((c(theta.new,beta.new)-c(theta.old,beta.old))**2)/sum(c(beta.old,theta.old))**2)
    theta.old<-theta.new
    beta.old<-beta.new
    p.old<-rep(theta.old[1], n)
    if (!is.null(x.0)){
      p.old<-rowSums(cbind(rep(theta.new[1], n),sapply(1:ncol(x.0), function(x) theta.new[x+1]*x.0[,x])))
    }
    if (is.finite(reldiff)==TRUE){
      if (reldiff<eps){
        conv<-TRUE
        break
      }
    }
  }
  list(conv=conv, reldiff=reldiff, it=it, mean=beta.new, theta.new=theta.new, fittedmean=mean.fit, p.old=p.old)
}