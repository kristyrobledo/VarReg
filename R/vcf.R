#' vcf dataset.
#'
#' A dataset containing 100 observations of mean velocity of circumferential fibre shortening (vcf), made by long axis and short axis echocardiography.
#'
#' @format A data frame with 133 rows and 3 variables:
#' \describe{
#'   \item{pid}{patient identifier}
#'   \item{vcflong}{vcf measurement from long axis}
#'   \item{vcfshort}{vcf measurement from short axis}
#'   ...
#' }
#' @source Data from Bland JM, Altman DG. (1986) Statistical methods for assessing agreement between two methods of clinical measurement. Lancet i, 307-310. (Supplied by Paul D'Arbela)


#'@examples
#'library(VarReg)
#'data(vcf)
#'attach(vcf)
#'plot(rowMeans(vcf[-1]),vcf$vcflong-vcf$vcfshort)
"vcf"