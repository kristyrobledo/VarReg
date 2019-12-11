#' SE calculations for mean and variance regression models
#'
#' \code{seVarReg} calculates SE for an object of class \code{VarReg}. If the result is not on a
#' boundary, the Fishers Information matrix SE are given. The bootstrapped 95\% CI can also be
#' calculated. Designed to be called by the plot function \code{plotVarReg}, rather than run by a user.
#' @param x Object of class \code{VarReg} to determin the SE (eg. result from \code{\link{semiVarReg}}).
#' @param boot Logical to indicate if bootstrapped CI should be calculated. Default is \code{FALSE}.
#' @param bootreps Number of bootstraps to be performed if \code{boot=TRUE}. Default is \code{1000}.
#' @param vector.mean Vector of \code{x} values for which the SE of the mean is to be calculated.
#' Default is the \code{x} covariate from the model.
#' @param vector.variance Vector of \code{x} values for which the SE of the variance is to be calculated.
#'  Default is the actual \code{x} covariate from the model.
#' @param control List of control parameters for the bootstrapped models.
#' See \code{\link{VarReg.control}}.
#' @param ... arguments to be used to form the default control argument if it is not supplied
#' directly
#' @return The result is a list of results. This includes:
#' \itemize{
#' \item \code{mean.est}: dataframe of overall results from the mean model, including parameter estimates
#' from the model, SEs from information matrix (if \code{boundary=FALSE}) and if specified, the SE
#' from bootstrapping with the bootstrapped 95\% CI.
#' \item \code{variance.est}: dataframe of overall results from the variance model, including parameter
#'  estimates from the model, SEs from information matrix (if \code{boundary=FALSE}) and if specified,
#'  the SE from bootstrapping with the bootstrapped 95\% CI.
#' \item \code{mean.im}: dataframe of the expected information matrices for the mean (as appropriate)
#' \item \code{variance.im}: dataframe of the expected information matrices for the variance
#' (as appropriate)
#' \item \code{mean.outputs}: dataframe with complete output for mean graphics. Includes the
#' \code{vector.mean} as input, and the mean vector (\code{mean.mean}) and the SE vector
#' \code{mean.se.im}, and bootstrapping outputs as appropriate.
#' \item \code{variance.outputs}: dataframe with complete output for variance graphics. Includes the
#' \code{vector.variance} as input, and the mean vector (\code{var.mean}) and the SE vector
#' \code{var.se.im}, and bootstrapping outputs as appropriate.
#' }
#'@seealso
#'\code{\link{semiVarReg}}, \code{\link{VarReg.control}}
#'
#'@examples
#'data(mcycle)
#'##Fit model with range as a covariate in the mean and the variance model
#'semimodel<-semiVarReg(mcycle$accel, mcycle$times, meanmodel="semi", varmodel="linear",
#'knots.m=4, maxit=10000)
#'##Calculate SE
#'se1<-seVarReg(semimodel, boot=FALSE)
#'##not run: with bootstrapping
#'##se2<-seVarReg(semimodel, boot=TRUE, bootreps=10)
#'##not run: calculate mean and SE for a given sequence
#'##test.seq<-seq(min(mcycle$times), max(mcycle$times),
#'##by=((max(mcycle$times)-min(mcycle$times))/999))
#'##se2<-seVarReg(semimodel, boot=TRUE, bootreps=10, vector.mean=test.seq)
#'@export

seVarReg<-function(x, boot=FALSE, bootreps=1000, vector.mean=x$data[,2],vector.variance=x$data[,2], control=list(...), ...){
  boot_est<-NULL
  n<-length(x$data[,1])
  control<-do.call(VarReg.control, control)
  if (x$boundary==TRUE){
    print("Note that as boundary==TRUE, SE must be obtained by bootstrapping")
    inv.beta<-NULL
    inv.alpha<-NULL

    if (x$meanmodel=="zero"){
      mean.mean<-rep(0, times=length(vector.mean))
      mean.d<-cbind(vector.mean, mean.mean, mean.se.im=NA)
    } else if (x$meanmodel=="constant"){
      pred.mean<-colSums(t(data.frame(rep(1, nrow(x$data)), x$data[,x$mean.ind]))*x$mean)
      mean.mean<-rep(pred.mean[1], times=length(vector.mean))
      mean.d<-cbind(vector.mean, mean.mean, mean.se.im=NA)
    } else if (x$meanmodel=="linear"){
      d<-cbind(rep(1, times=length(vector.mean)), vector.mean)
      mean.mean<-colSums(t(d)*x$mean)
      mean.d<-cbind(vector.mean, mean.mean, mean.se.im=NA)
    }else if (x$meanmodel=="semi"){
      bmean<-splines::bs(x$data[,2], degree=x$degree, df=x$degree+x$knots.m)
      k <- attr(bmean, "knots")
      d<-cbind(rep(1, times=length(vector.mean)), splines::bs(vector.mean, degree=x$degree, knots=k))
      mean.mean<-colSums(t(d)*x$mean)
      mean.d<-cbind(vector.mean, mean.mean, mean.se.im=NA)
    }

    if (x$varmodel=="constant"){
      pred.var<-colSums(t(data.frame(rep(1, nrow(x$data)), x$data[,x$var.ind]))*x$variance)
      var.mean<-rep(pred.var[1], times=length(vector.variance))
      var.d<-cbind(vector.variance, var.mean, var.se.im=NA)
    } else if (x$varmodel=="linear"){
      d2<-cbind(rep(1, times=length(vector.variance)), vector.variance)
      var.mean<-colSums(t(d2)*x$variance)
      var.d<-cbind(vector.variance, var.mean, var.se.im=NA)
    }else if (x$varmodel=="semi"){
      bvar<-splines::bs(x$data[,2], degree=x$degree, df=x$degree+x$knots.v)
      k <- attr(bvar, "knots")
      d2<-cbind(rep(1, times=length(vector.variance)), splines::bs(vector.variance, degree=x$degree, knots=k))
      var.mean<-colSums(t(d2)*x$variance)
      var.d<-cbind(vector.variance, var.mean, var.se.im=NA)
    }


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
      var.mean<-rep(pred.var[1], times=length(vector.variance))
      var.se.im<-rep(sqrt(inv.alpha), times=length(vector.variance) )
      var.d<-cbind(vector.variance, var.mean, var.se.im)
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
      if (x$varmodel=="linear"){
        d2<-cbind(rep(1, times=length(vector.variance)), vector.variance)
      }else if (x$varmodel=="semi"){
        bvar<-splines::bs(x$data[,2], degree=x$degree, df=x$degree+x$knots.v)
        k <- attr(bvar, "knots")
        d2<-cbind(rep(1, times=length(vector.variance)), splines::bs(vector.variance, degree=x$degree, knots=k))
      }
      var.mean<-colSums(t(d2)*x$variance)
      var.var.im<-vector(length=length(vector.variance))
      for (i in 1:ncol(d2)){
        for (j in 1:ncol(d2)){
          var.var.im<-var.var.im+d2[,i]*d2[,j]*inv.alpha[i,j]
        }
      }
      var.se.im<-sqrt(var.var.im)
      var.d<-cbind(vector.variance, var.mean, var.se.im)
    }


    if (x$meanmodel=="zero"){
      inv.beta=NULL
      mean.mean<-rep(0, length(vector.mean))
      mean.d<-cbind(vector.mean, mean.mean, mean.se.im=NA)
    }else if (x$meanmodel=="constant"){
      pred.mean<-colSums(t(data.frame(rep(1, nrow(x$data)), x$data[,x$mean.ind]))*x$mean)
      exp.i.mean<- sum(1/(pred.var))
      inv.beta<-tryCatch({
        solve(exp.i.mean, tol=1e-28)
      },error=function(cond){
        message("Failed to calculate the Inverse expected I-matrix for the mean")
        return(matrix(NA, nrow=(length(x$mean.ind)+1), ncol=(length(x$mean.ind)+1)))
      })
      mean.mean<-rep(pred.mean[1], times=length(vector.mean))
      mean.se.im<-rep(sqrt(inv.beta), times=length(vector.mean) )
      mean.d<-cbind(vector.mean, mean.mean, mean.se.im)
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

      if (x$meanmodel=="linear"){
        d<-cbind(rep(1, times=length(vector.mean)), vector.mean)
      }else if (x$meanmodel=="semi"){
        bmean<-splines::bs(x$data[,2], degree=x$degree, df=x$degree+x$knots.m)
        k <- attr(bmean, "knots")
        d<-cbind(rep(1, times=length(vector.mean)), splines::bs(vector.mean, degree=x$degree, knots=k))
      }
      mean.mean<-colSums(t(d)*x$mean)
      mean.var.im<-vector(length=length(vector.mean))
      for (i in 1:ncol(d)){
        for (j in 1:ncol(d)){
          mean.var.im<-mean.var.im+d[,i]*d[,j]*inv.beta[i,j]
        }
      }
      mean.se.im<-sqrt(mean.var.im)
      mean.d<-cbind(vector.mean, mean.mean, mean.se.im)
    }
  }
  if (boot==TRUE){
    m.reps<-matrix(NA,bootreps,length(x$mean),dimnames=list(seq(1,bootreps, by=1),names(x$mean)))
    v.reps<-matrix(NA, bootreps,length(x$variance),dimnames=list(seq(1,bootreps, by=1),names(x$variance)))
    bases_mean<-matrix(NA, nrow=bootreps, ncol=length(vector.mean))
    bases_var<-matrix(NA, nrow=bootreps, ncol=length(vector.variance))
    for (b in 1:bootreps){
      drep<-sample(1:n, replace=TRUE)
      Yrep<-x$data[drep,1]
      xrep<-x$data[drep,2]
      if (is.null(x$cens.ind)==FALSE){
        ccrep<-x$data[drep,x$cens.ind]
      }else if (is.null(x$cens.ind)==TRUE){
        ccrep<-NULL
        }
      log <- utils::capture.output({
        bootmodel<-semiVarReg(y=Yrep,x=xrep, cens.ind=ccrep, meanmodel=x$meanmodel, varmodel=x$varmodel, knots.m=x$knots.m, knots.v=x$knots.v, degree=x$degree, mono.var=x$mono.var, control=control)
      })
      m.reps[b,]<-bootmodel$mean
      v.reps[b,]<-bootmodel$variance
      if (x$meanmodel=="constant"){
        bases_mean[b,]<-rep(bootmodel$mean, times=length(vector.mean))
      }else if (x$meanmodel=="linear"){
        bases_mean[b,]<-colSums(t(d)*bootmodel$mean)
      }else if (x$meanmodel=="semi"){
        b1<-splines::bs(xrep, degree=x$degree, df=(x$degree+x$knots.m))
        b2<-splines::bs(vector.mean, degree=x$degree, knots=attr(b1, "knots"), Boundary.knots = c(min(attr(b1, "Boundary.knots")[1], vector.mean), max(attr(b1, "Boundary.knots")[2], vector.mean)))
        bases_mean[b,]<-colSums(t(cbind(rep(1,length(vector.mean)), b2))*bootmodel$mean)
      }
      if (x$varmodel=="constant"){
        bases_var[b,]<-rep(bootmodel$variance, times=length(vector.variance))
      }else if (x$varmodel=="linear"){
        bases_var[b,]<-colSums(t(d2)*bootmodel$variance)
      }else if (x$varmodel=="semi"){
        b1v<-splines::bs(xrep, degree=x$degree, df=(x$degree+x$knots.v))
        b2v<-splines::bs(vector.variance, degree=x$degree, knots=attr(b1v, "knots"), Boundary.knots = c(min(attr(b1v, "Boundary.knots")[1], vector.variance), max(attr(b1v, "Boundary.knots")[2], vector.variance)))
        bases_var[b,]<-colSums(t(cbind(rep(1,length(vector.variance)), b2v))*bootmodel$variance)
      }
    }
    if (x$meanmodel!="zero" ){
  mean.d<-cbind(mean.d, t(apply(bases_mean , 2 ,stats::quantile, probs=c(0.50,0.025, 0.975),na.rm=TRUE)))
  colnames(mean.d)[4:6]<-c("boot.median", "boot.lci", "boot.uci")
    }
    var.d<-cbind(var.d, t(apply(bases_var , 2 ,stats::quantile, probs=c(0.50,0.025, 0.975),na.rm=TRUE)))
    colnames(var.d)[4:6]<-c("boot.median", "boot.lci", "boot.uci")
  }

  if (x$meanmodel=="zero"){
      m.estimates<-NULL
  }else {
      m.estimates<-matrix(NA, nrow=length(x$mean), ncol=5, dimnames=list(names(x$mean), c("Estimate", "se.im", "se.boot", "lci.boot", "uci.boot")))
      m.estimates[,1]<-x$mean
      if (is.null(inv.beta)==FALSE){
        m.estimates[,2]<-diag(inv.beta)
      }
      if (boot==TRUE && x$meanmodel!="semi"){
        m.estimates[,3]<-apply(m.reps, 2, function(x) stats::sd(x, na.rm=TRUE))
        m.estimates[,4]<-apply(m.reps, 2, function(x) stats::quantile(x, prob=0.025, na.rm=TRUE))
        m.estimates[,5]<-apply(m.reps, 2, function(x) stats::quantile(x, prob=0.975, na.rm=TRUE))
    }
  }

  v.estimates<-matrix(NA, nrow=length(x$variance), ncol=5, dimnames=list(names(x$variance), c("Estimate", "se.im", "se.boot", "lci.boot", "uci.boot")))
  v.estimates[,1]<-x$variance
  if (is.null(inv.alpha)==FALSE){
    v.estimates[,2]<-diag(inv.alpha)
  }
  if (boot==TRUE && x$varmodel!="semi"){
    v.estimates[,3]<-apply(v.reps, 2, function(x) stats::sd(x, na.rm=TRUE))
    v.estimates[,4]<-apply(v.reps, 2, function(x) stats::quantile(x, prob=0.025, na.rm=TRUE))
    v.estimates[,5]<-apply(v.reps, 2, function(x) stats::quantile(x, prob=0.975, na.rm=TRUE))
  }


  list(mean.est=m.estimates, variance.est=v.estimates, mean.im=inv.beta, variance.im=inv.alpha, mean.outputs=mean.d, variance.outputs=var.d)
}



