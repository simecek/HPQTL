#' Half-inverse matrix operation
#'  
#' @param W matrix
#'
#' @references Abney M, Ober C, McPeek MS (2002). "Quantitative trait homozygosity
#' and association mapping and empirical genome-wide significance in
#' large complex pedigrees: Fasting serum insulin levels in the
#' Hutterites." American Journal of Human Genetics 70: 920-934.
#'
#' @keywords manip
#' 
#' @value list
#'
#' @author copied from QTLRel
#'
#' @examples
#' cross <- sim.cross.geno(250, nmar=10)
#' cross$pheno <- sim.cross.pheno(0.5, cross)
#' heritability(cross)

half.inv <- function(W, symmetric=TRUE,inverse=TRUE){
  eW <- eigen(W, symmetric=symmetric)
  d <- eW$values
  if (min(d) <0  && abs(min(d))>sqrt(.Machine$double.eps))
    stop("'W' is not positive definite")
  else d[d<=0]<- ifelse(inverse, Inf, 0)
  A <- diag(d^ifelse(inverse, -0.5, 0.5)) %*% t(eW$vector)
  A # t(A)%*%A = W^{-1}
}