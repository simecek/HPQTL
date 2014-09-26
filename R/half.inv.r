#' Half-inverse matrix operation
#'
#' Return square root of an inverse matrix. 
#' 
#' @param W positive definite matrix
#'
#' @details Note \code{W} must be strictly positive definite. 
#' For a matrix with negative eigen values the error is thrown. However, the problem
#' might arise for matrices with eigen values close to zero. 
#'
#' @keywords manip
#' 
#' @return Matrix \code{A} such that \code{t(A)\%*\%A} is the inverse of \code{W}. 
#'
#' @author Riyan Cheng, copied from QTLRel
#'
#' @examples
#' W <- matrix(c(1,1/2,1/2,1),2,2)
#' A <- HPQTL:::half.inv(W)
#' max(abs(t(A)%*%A%*%W) - diag(2))

half.inv <- function(W, symmetric=TRUE, inverse=TRUE) {
  eW <- eigen(W, symmetric=symmetric)
  d <- eW$values
  if (min(d)<0  && abs(min(d))>sqrt(.Machine$double.eps))
    stop("'W' is not positive definite")
  else d[d<=0] <- ifelse(inverse, Inf, 0)
  A <- diag(d^ifelse(inverse, -0.5, 0.5)) %*% t(eW$vector)
  A # t(A)%*%A = W^{-1}
}