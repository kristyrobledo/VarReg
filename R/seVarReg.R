#' SE for mean and variance regression
#'
#' \code{seVarReg} calculates SE for an object from a linVarReg. If the result is not on a boundary, the Fishers Information matrix SE are given. The bootstrapped 95% CI can also be calculated.
#' @param x object of class VarReg to determin the SE (eg. model from semiVarReg)
#' @param boot TRUE or FALSE indicating if bootstrapped CI should be calculated. Default is FALSE.
#' @param bootreps Number of bootstraps to be performed if $boot=TRUE$. Default is 1000.
#' @return Table of results.
#'@examples
#'data(lidar)
#'lid<-data.frame(lidar$logratio, lidar$range)
#'##Fit model with range as a covariate in the mean and the variance model
#'semimodel<-semiVarReg(lidar$logratio, lidar$range, meanmodel="semi", varmodel="semi", knots.m=4, knots.v=2)
#'##Calculate SE
#'se<-seVarReg(semimodel, boot=TRUE, bootreps=10)
#'@export

seVarReg<-function(x, boot=FALSE, bootreps=1000, eps=1e-6, maxit=1000){
  n<-length(x$data[,1])
  if (x$boundary==TRUE){
    print("As boundary==TRUE, SE must be obtained by bootstrapping")
  }else if (x$boundary==FALSE){
    pred.var<-colSums(t(data.frame(rep(1, nrow(x$data)), x$data[,x$var.ind]))*x$variance)
    if (x$varmodel=="constant"){
      exp.i.var<- 1/2*sum(1/(pred.var**2))
      inv.alpha<-tryCatch({
        solve(exp.i.var, tol=1e-28)
      },error=function(cond){
        message("Failed to calculate the Inverse expected I-matrix for the variance")
        return(matrix(NA, nrow=(length(x$var.ind)+1), ncol=(length(x$var.ind)+1)))
      })
    }else {
      #matrix for variance
      exp.i.var<-matrix(NA, nrow=(length(x$var.ind)+1), ncol=(length(x$var.ind)+1))

      for (j in 1:(length(x$var.ind)+1)){
        for (k in 1:(length(x$var.ind)+1)){
          if (j==1){
            if (k==1){
              ##intercept
              exp.i.var[1,1]<- 1/2*sum(1/(pred.var**2))
            }else{
              ##first column
              exp.i.var[k,1]<- 1/2*sum(x$data[,x$var.ind[k-1]]/(pred.var**2))
            }
          }else if (k==1){
            ##first row
            exp.i.var[1,j]<- 1/2*sum(x$data[,x$var.ind[j-1]]/(pred.var**2))
          }else{
            #remains of the matrix
            exp.i.var[k,j]<-1/2*sum(x$data[,x$var.ind[k-1]]*(x$data[,x$var.ind[j-1]])/(pred.var**2))
          }
        }
      }

      inv.alpha<-tryCatch({
        solve(exp.i.var, tol=1e-28)
      },error=function(cond){
        message("Failed to calculate the Inverse expected I-matrix for the variance")
        return(matrix(NA, nrow=(length(x$var.ind)+1), ncol=(length(x$var.ind)+1)))
      })
    }
    if (x$meanmodel=="zero"){
      pred.mean<-rep(0, nrow(x$data))
      inv.beta=NULL
    }else if (x$meanmodel=="constant"){
      pred.mean<-colSums(t(data.frame(rep(1, nrow(x$data)), x$data[,x$mean.ind]))*x$mean)
      exp.i.mean<- sum(1/(pred.var))
      inv.beta<-tryCatch({
        solve(exp.i.mean, tol=1e-28)
      },error=function(cond){
        message("Failed to calculate the Inverse expected I-matrix for the mean")
        return(matrix(NA, nrow=(length(x$mean.ind)+1), ncol=(length(x$mean.ind)+1)))
      })
    }else {
      pred.mean<-colSums(t(data.frame(rep(1, nrow(x$data)), x$data[,x$mean.ind]))*x$mean)
      exp.i.mean<-matrix(NA, nrow=(length(x$mean.ind)+1), ncol=(length(x$mean.ind)+1))
      for (j in 1:(length(x$mean.ind)+1)){
        # print(c("J", j))
        for (k in 1:(length(x$mean.ind)+1)){
          # print(c("K", k))
          if (j==1){
            if (k==1){
              exp.i.mean[1,1]<- sum(1/(pred.var))
            }else{
              exp.i.mean[k,1]<- sum(x$data[,x$mean.ind[k-1]]/pred.var)
            }
          }else if (k==1){
            exp.i.mean[1,j]<- sum(x$data[,x$mean.ind[j-1]]/pred.var)
          }else{
            exp.i.mean[k,j]<- sum((x$data[,x$mean.ind[k-1]]*x$data[,x$mean.ind[j-1]])/pred.var)
          }
        }
      }
      inv.beta<-tryCatch({
        solve(exp.i.mean, tol=1e-28)
      },error=function(cond){
        message("Failed to calculate the Inverse expected I-matrix for the mean")
        return(matrix(NA, nrow=(length(x$mean.ind)+1), ncol=(length(x$mean.ind)+1)))
      })
    }
  }
  if (boot==TRUE){
    m.reps<-matrix(NA,bootreps,length(x$mean),dimnames=list(seq(1,bootreps, by=1),names(x$mean)))
    v.reps<-matrix(NA, bootreps,length(x$variance),dimnames=list(seq(1,bootreps, by=1),names(x$variance)))
    for (b in 1:bootreps){
      drep<-sample(1:n, replace=TRUE)
      Yrep<-x$data[drep,1]
      if (x$meanmodel=="linear"){
      xrep<-x$data[drep,2]
      } else if (x$meanmodel=="semi"){
      xrep<-x$data[drep,2]
      }
      log <- capture.output({
        bootmodel<-semiVarReg(y=Yrep,x=xrep, meanmodel=x$meanmodel, varmodel=x$varmodel, knots.m=x$knots.m, knots.v=x$knots.v, degree=x$degree, mono.var=x$mono.var, eps=eps, maxit=maxit)
      })
      m.reps[b,]<-bootmodel$mean
      v.reps[b,]<-bootmodel$variance
    }
  }
  if (x$meanmodel=="zero"){
      m.estimates<-NULL
  }else{
      m.estimates<-matrix(NA, nrow=length(x$mean), ncol=5, dimnames=list(names(x$mean), c("Estimate", "se.im", "se.boot", "lci.boot", "uci.boot")))
      m.estimates[,1]<-x$mean
      m.estimates[,2]<-diag(inv.beta)
      if (boot==TRUE){
        m.estimates[,3]<-apply(m.reps, 2, function(x) sd(x, na.rm=TRUE))
        m.estimates[,4]<-apply(m.reps, 2, function(x) quantile(x, prob=0.025, na.rm=TRUE))
        m.estimates[,5]<-apply(m.reps, 2, function(x) quantile(x, prob=0.975, na.rm=TRUE))
    }
  }

  v.estimates<-matrix(NA, nrow=length(x$variance), ncol=5, dimnames=list(names(x$variance), c("Estimate", "se.im", "se.boot", "lci.boot", "uci.boot")))
  v.estimates[,1]<-x$variance
  v.estimates[,2]<-diag(inv.alpha)
  if (boot==TRUE){
    v.estimates[,3]<-apply(v.reps, 2, function(x) sd(x, na.rm=TRUE))
    v.estimates[,4]<-apply(v.reps, 2, function(x) quantile(x, prob=0.025, na.rm=TRUE))
    v.estimates[,5]<-apply(v.reps, 2, function(x) quantile(x, prob=0.975, na.rm=TRUE))
  }
  list(mean.est=m.estimates, variance.est=v.estimates, mean.im=inv.beta, variance.im=inv.alpha)
}



