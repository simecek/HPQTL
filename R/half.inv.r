#' Half-inverse matrix operation
#'  
#' @param W matrix
#'
#' @keywords manip
#' 
#' @return matrix
#'
#' @author copied from QTLRel
#'
#' @examples
#' data(fake.f2, package="qtl")
#' fake.f2 <- calc.genoprob(fake.f2)
#' 
#' geno <- extract.geno(fake.f2)
#' G <- gensim.matrix(geno)
#' A <- HPQTL:::half.inv(G)

half.inv <- function(W, symmetric=TRUE,inverse=TRUE){
  eW <- eigen(W, symmetric=symmetric)
  d <- eW$values
  if (min(d) <0  && abs(min(d))>sqrt(.Machine$double.eps))
    stop("'W' is not positive definite")
  else d[d<=0]<- ifelse(inverse, Inf, 0)
  A <- diag(d^ifelse(inverse, -0.5, 0.5)) %*% t(eW$vector)
  A # t(A)%*%A = W^{-1}
}