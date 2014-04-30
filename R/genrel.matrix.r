#' Genetic Relationship Matrix
#' 
#' 
#' @param geno \code{genotype.probs} object or \code{cross} object
#' @param method how to calculate similarity matrix ('allele-2f-additive', 'allele-multif-additive')    
#' @exclude.chrs
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
#' genrel.matrix(geno, method = "allele-multif-additive")

genrel.matrix <- function(geno, method=c('allele-2f-additive', 'allele-multif-additive'), skip.chrs=c()) {
  
  method = match.arg(method)
  
  # if not genotype.probs then extract genotype probs
  if ("cross" %in% class(geno)) geno <- extract.geno(geno)
   
  if (method == 'allele-2f-additive') {
    score.matrix <- matrix(c(1,1/2,0,1/2,1,1/2,0,1/2,1),3,3)
    
    K <- matrix(0, nrow(geno$probs[[1]]), nrow(geno$probs[[1]]))
    for (i in setdiff(1:length(geno$probs),skip.chrs))
      if (length(geno$calls[[i]])==3) { # AA,AB,BB calls
        for (j in 1:dim(geno$probs[[i]])[3]) # for each snp get kinship
          K <- K + geno$probs[[i]][,,j] %*% score.matrix %*% t(geno$probs[[i]][,,j])
      } else { # AA,AB calls
        for (j in 1:dim(geno$probs[[i]])[3]) # for each snp get kinship
          K <- K + geno$probs[[i]][,,j] %*% score.matrix[1:2,1:2] %*% t(geno$probs[[i]][,,j])
      }    
      
    # normalization to diag(K) == 1
    norm.factor <- sqrt(diag(K))
    K <- diag(1/norm.factor) %*% K %*%  diag(1/norm.factor)
    
    return(K)
  }
  
  if (method == 'allele-multif-additive') {
    
    K <- matrix(0, nrow(geno$probs[[1]]), nrow(geno$probs[[1]]))
    for (i in setdiff(1:length(geno$probs),skip.chrs)) {
      score.matrix <- diag(length(geno$calls[[i]]))
      for (j in 1:dim(geno$probs[[i]])[3]) # for each snp get kinship
        K <- K + geno$probs[[i]][,,j] %*% score.matrix %*% t(geno$probs[[i]][,,j])
    }    
    
    # normalization to diag(K) == 1
    norm.factor <- sqrt(diag(K))
    K <- diag(1/norm.factor) %*% K %*%  diag(1/norm.factor)
    
    return(K)
  }
  
  stop("The method is not yet implemented")  
  
}