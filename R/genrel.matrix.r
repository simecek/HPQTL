#' Genetic Relationship Matrix
#' 
#' 
#' @param geno \code{cross} object or a genotype matrix (\code{nind} rows, \code{totalmar} columns)
#' @param method how to calculate similarity matrix ('additive', 'additive-diag1', 'kinship')
#' @param expected.means a vector of expected marker means, if \code{NULL} means are estimated from data    
#'
#' @return square matrix (\code{nind} rows, \code{nind} columns)
#' 
#' @keywords manip
#'
#' @export
#' 
#' @examples
#' cross <- sim.cross.geno(250, nmar=10)
#' geno <- extract.geno(cross)
#' genrel.matrix(geno)
#' genrel.matrix(geno, method = "kinship")
#' genrel.matrix(geno, expected.means = c(rep(2,190), rep(1.5,10)))

genrel.matrix <- function(geno, method="additive", expected.means=NULL) {
  
  # if not matrix then extract genotype
  if ("cross" %in% class(geno)) geno <- extract.geno(geno)
  
  # if means not give, estimate them
  if (is.null(expected.means)) expected.means <- colMeans(geno) 
  
  if (method == "kinship") {
    K <- diag(rep(1, nrow(geno)))
    for (i in 1:(nrow(geno)-1))
      for (j in (i+1):nrow(geno))
        K[i,j] <- K[j,i] <- 1 - mean(abs(geno[i,]-geno[j,]))/2
    return(K)
  }
  
  if (method == "additive" || method == "additive-diag1") {
    
    G <- (geno - expected.means) %*% t(geno - expected.means)
    
    # normalization
    if (method == "additive") {
      G <- G / mean(diag(G))
    } else {
      norm.factor <- sqrt(diag(G))
      G <- diag(1/norm.factor) %*% G %*%  diag(1/norm.factor)
    }
    
    return(G)
  }
  
}