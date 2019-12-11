#' Calculation of information criterion
#'
#' \code{criterion} calculates various information criterion for the algorithms in this package
#' @param n Number of observations
#' @param loglik Loglikelihood from model
#' @param param Number of parameters fit in model
#' @return A list of the four IC
#' \itemize{
#'  \item\code{aic.c}: Akaike information criterion corrected for small samples
#'  \item\code{aic}: Akaike information criterion
#'  \item\code{bic}: Bayesian information criterion
#'  \item\code{hqc}: Hannan-Quinn information criterion
#'  }
#' @export

criterion<-function(n, loglik, param){
  aicc<-(2*(param)*n)/(n-(param)-1)-2*loglik
  aic<-2*(param)-2*loglik
  bic<-log(n)*(param)-2*loglik
  hqc<-log(log(n))*(param)*2-2*loglik
  list(aicc=aicc, aic=aic, bic=bic, hqc=hqc)
}