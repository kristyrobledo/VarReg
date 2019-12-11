#' The EM loop for the LSS model
#'
#' \code{loop_lss} is the EM loop function for the LSS model to be utilised by various other higher level functions
#' @param alldat Dataframe containing all the data for the models. Outcome in the first column.
#' @param xiold Vector of initial location parameter estimates to be fit in the location model.
#' @param omega2old Vector of initial scale2 parameter estimates to be fit in the scale2 model.
#' @param nuold Vector of initial nu parameter estimates to be fit in the nu model.
#' @param mean.ind Vector containing the column numbers of the data in 'alldat' to be fit as covariates in the location model.
#' @param var.ind Vector containing the column numbers of the data in 'alldat' to be fit as covariates in the scale2 model. FALSE indicates a constant variance model.
#' @param nu.ind Vector containing the column numbers of the data in 'alldat' to be fit as covariates in the nu model. NULL indicates constant model.
#' @param para.space Parameter space to search for variance parameter estimates. "positive" means only search positive parameter space, "negative" means search only negative parameter space and "all" means search all.
#' @param maxit Number of maximum iterations for the main EM algorithm.
#' @param eps Very small number for the convergence criteria.
#' @param int.maxit Number of maximum iterations for the internal EM algorithm for the location and scale.
#' @param print.it Logical to indicate if the estimates for each iteration should be printed.
#' @return A list of the results from the algorithm, including conv, reldiff, it, mean, xi.new, omega2.new, nu.new, fitted.xi
#' \itemize{
#' \item\code{conv}: Logical argument indicating if convergence occurred
#' \item\code{it}: Total iterations performed of the EM algorithm
#'  \item\code{reldiff}: the positive convergence tolerance that occured at the final iteration
#'  \item\code{xinew}: Vector of location parameter estimates
#'  \item\code{omega2new}: Vector of scale squared parameter estimates
#'  \item\code{nunew}: Vector of shape parameter estimates
#'  \item\code{fitted.xi}: Vector of fitted location estimates
#'  }
#'@export

loop_lss<-function(alldat,xiold,omega2old,nuold,mean.ind, var.ind, nu.ind, para.space,maxit,eps,int.maxit, print.it){
  conv<-FALSE
  Y<-alldat[,1]
  n<-length(Y)
  cc<-c(rep(0,n),rep(-1,n))
  thetaold<-c(xiold,omega2old,nuold)
  for (a in 1:maxit) {
    mean<-unname(colSums(t(alldat[1:n,mean.ind])*xiold))
    variance<-unname(colSums(t(cbind(rep(1,n),alldat[1:n,var.ind]))*omega2old))
    r<- ((Y-mean)/sqrt(variance))
    r[round(r, digits=5)==0]<-0.0000001
    yy<-r>=0
    if (length(nuold)==1){
      nunew<-unname(stats::glm(yy~-1+abs(r),family=stats::binomial(link="probit"),
                          stats::glm.control(epsilon = eps, maxit = maxit))$coef)    # no intercept
      nu<-rep(nunew, n)
    }else{
      ##regression model in scale as well
      reg<-data.frame(apply(cbind(rep(1,n), alldat[,nu.ind]), 2, function(x) x*abs(r)))
      gg<-stats::glm(yy~-1+. , data=reg, family=stats::binomial(link="probit"), epsilon = eps, maxit = maxit)
      nunew<-unname(gg$coef)
      nu<-unname((gg$linear.predictors)/abs(r))
    }
    if (var.ind[1]==FALSE){
      dat2<-rbind(alldat, apply(alldat, 2, function(x) x*nu))
      dat3<-cbind(dat2, cc)
    }else{
      dat2<-rbind(alldat[-var.ind], apply(alldat[-var.ind], 2, function(x) x*nu))
      dat3<-cbind(dat2,rbind(as.data.frame(alldat[,var.ind]), as.data.frame(alldat[,var.ind])), cc)
    }

    cc.ind<-which(colnames(dat3)%in%"cc")
    if (print.it==FALSE){
      log <- utils::capture.output({
      output<-censlinVarReg(dat=dat3,mean.ind=mean.ind,var.ind=var.ind,cens.ind=cc.ind, mean.intercept=FALSE, para.space=para.space, mean.init = xiold, var.init=omega2old, maxit=int.maxit)
       })
    }else if (print.it==TRUE){
      print(paste0("Iteration: ", a))
      output<-censlinVarReg(dat=dat3,mean.ind=mean.ind,var.ind=var.ind,cens.ind=cc.ind, mean.intercept=FALSE, para.space=para.space, mean.init = xiold, var.init=omega2old, maxit=int.maxit)
      print(c("Location: ", paste(output$mean)))
      print(c("Scale2: ", paste(output$variance)))
      print(c("Shape: ", paste(nunew)))
    }

    if (is.na(output$reldiff)==TRUE){
      message("Location and scale model has not converged. Please review initial estimates.")
      reldiff<-NA
      break
    }
    xinew<-output$mean
    omega2new<-output$variance
    thetanew<-c(xinew, omega2new, nunew)
    pdiff<-sqrt(sum((thetanew-thetaold)**2)/sum(thetaold**2))
    if (pdiff<eps) {
      conv<-TRUE
      break
    }
    xiold<-xinew
    omega2old<-omega2new
    nuold<-nunew
    thetaold<-thetanew
  }
  mean<-unname(colSums(t(alldat[1:n,mean.ind])*xiold))
  list(conv=conv, reldiff=pdiff, it=a, xinew=xinew, omega2new=omega2new, nunew=nunew, fitted.xi=mean)
}