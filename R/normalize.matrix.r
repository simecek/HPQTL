#' Normalize Matrix
#'
#' Normalize genetic similarity matrix.
#' 
#' @param W symmetric, positive semidefinite matrix
#' @param method "sample-variance", "diagonal" or "diagonal-average"
#'
#' @details To calculate heritability, particularly for different genetic similarity estimators,
#' the genetic similarity matrix needs to be properly normalized. The default way is to force sample
#' variance of animal effect to be unit. 
#'
#' @keywords manip
#' 
#' @return Normalized matrix
#'
#' @examples
#' W <- matrix(c(1,1/2,1/2,1),2,2)
#' normalize.matrix(W, method="sample-variance")
#' normalize.matrix(W, method="diagonal")

normalize.matrix <- function(W, method=c("sample-variance", "diagonal", "diagonal-average", "none")) {

  # check parameters
  method <- match.arg(method)
  stopifnot(class(W)=="matrix" | class(W)=="list")
  
  if (class(W)=="list")
    lapply(W, normalize.matrix, method=method)
  
  stopifnot(is.numeric(W) & is.finite(W))
  
  if (method == "sample-variance") {
     n <- nrow(W)
     k <- 1/n * (sum(diag(W)) - 2*sum(W[upper.tri(W)])/(n-1))
     return(W/k)
  }
  
  # make all diagonal elements equal one
  if (method == "diagonal") {
     V <- diag(1/sqrt(diag(W)))
     return(V%*%W%*%V)
  }
  
  # make average diagonal element equals one
  if (method == "diagonal-average") {
     k <- mean(diag(W))
     return(W/k)
  }
  
  if (method == "none") {
    return(W)
  }  
}