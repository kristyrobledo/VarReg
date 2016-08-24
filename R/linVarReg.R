#' Linear mean and variance regression
#'
#' \code{linVarReg} performs multivariate mean and multivariate variance regression.
#' @param dat Dataframe containing outcome and covariate data. Outcome data must be in the first column. Covariates for mean and variance model in next columns.
#' @param var.ind Vector containing the column numbers of the data in 'dat' to be fit as covariates in the variance model. FALSE indicates constant variance option.
#' @param mean.ind Vector containing the column numbers of the data in 'dat' to be fit as covariates in the mean model. 0 indicates constant mean option. NULL indicates zero mean option.
#' @param para.space Parameter space to search for variance parameter estimates. "positive" means only search positive parameter space, "negative" means search only negative parameter space and "all" means search all.
#' @param maxit Number of maximum iterations for the EM algorithm, default 1000.
#' @param eps Very small number for the convergence criteria, default 1 times 10 power6
#' @return $linVarReg$ returns an object of class "VarReg" which inherits some components from the class "glm".
#'
#' A list of the results from the algorithm, including conv, reldiff, information criterion and mean and variance estimates.
#'@examples
#'data(lidar)
#'lid<-data.frame(lidar$logratio, lidar$range)
#'##Fit model with range as a covariate in the mean and the variance model
#'linmodel<-linVarReg(dat=lidar, var.ind=c(2),mean.ind=c(2), para.space="all")
#'##Give more iterations to converge
#'linmodel<-linVarReg(dat=lidar, var.ind=c(2),mean.ind=c(2), maxit=10000, para.space="all")
#'##Fit model with range as variate in variance model only - constant mean model.
#'conmean<-linVarReg(dat=lidar, var.ind=c(2),mean.ind=0, para.space="all")
#'@export


linVarReg<-function(dat, var.ind=c(2), mean.ind=c(2), para.space=c("all", "positive", "negative"), eps=1e-6, maxit=1000){
  para.space<-match.arg(para.space)
  loops<-list()
  ll<-vector()
   X<-dat[,1]
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
    loops[[1]]<-loop_em(meanmodel, theta.old=1, p.old=rep(1, n), x.0=NULL, X, maxit, eps)
    var<- rep(loops[[1]]$theta.new, n)
    ll<- -n/2*log(2*pi)-1/2*sum(log(var))-sum(((X-loops[[1]]$fitted)**2)/(2*(var)))
    if (loops[[1]]$conv==FALSE){
      writeLines(paste("Warning: Did not converge at Maxit=", maxit))
    }

  ##variance model option
  }else{
    ## set parameters
    minmax<-list()
    x<-as.matrix(dat[,var.ind])
    theta.old<-rep(1, 1+length(var.ind))
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
      loops[[i]]<-loop_em(meanmodel, theta.old, p.old, x.0, X, maxit, eps)

      ll[i]<- -n/2*log(2*pi)-1/2*sum(log(loops[[i]]$p.old))-sum(((X-loops[[i]]$fittedmean)**2)/(2*(loops[[i]]$p.old)))
      if (loops[[i]]$conv==FALSE){
        writeLines(paste("Warning: Did not converge at maxit=", maxit))
      }
    }
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

    names(alpha)<-c("Intercept", colnames(dat)[var.ind])
   if (is.null(mean.ind)==FALSE){
     names(mean)<-c("Intercept", colnames(dat)[mean.ind])
   }
  if (is.null(mean.ind[1])==TRUE){
    #print("test param")
    param<-length(totalx)+1
  }else if (mean.ind[1]==0){
    param<-1+length(totalx)+1
  }else{
    param<-length(mean.ind)+1+length(totalx)+1
  }

  if (sum(as.integer(loops[[max.ll]]$p.old < eps))>0){
    boundary=TRUE
  }else{
    boundary=FALSE
  }
  aicc<-(2*(param)*n)/(n-(param)-1)-2*ll.final
  aic<-2*(param)-2*ll.final
  bic<-log(n)*(param)-2*ll.final
  hqc<-log(log(n))*(param)-2*ll.final
  fit<- list(converged=conv.final,iterations=it.final,reldiff=reldiff.final, loglik=ll.final, boundary=boundary, aic.c=aicc, aic=aic,bic=bic,hqc=hqc,mean.ind=mean.ind,mean=mean,var.ind=var.ind, variance=alpha, data=dat)
  class(fit) <- c("VarReg")
  return(fit)
}
