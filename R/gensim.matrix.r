#' Genetic Similarity Matrix
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
#' gensim.matrix(geno)
#' gensim.matrix(geno, method = "allele-multif-additive")

gensim.matrix <- function(geno, method=c('allele-2f-additive', 'allele-multif-additive', 'allele-multif-cosine'), 
                          skip.chrs=c()) {
  
  method = match.arg(method)
  
  # if not genotype.probs then extract genotype probs
  if ("cross" %in% class(geno)) geno <- extract.geno(geno)

  # select all snps but those from chromosome skip.chrs
  if (length(skip.chrs)>0) {
    sel.snps <- which(!(geno$markers$chr %in% skip.chrs))
  } else {
    sel.snps <- 1:dim(geno$probs)[[3]] # all
  }
  
  if (method == 'allele-2f-additive') {
    score.matrix <- matrix(c(1,1/2,0,1/2,1,1/2,0,1/2,1),3,3)
    K <- matrix(0, nrow(geno$probs), nrow(geno$probs))
    
    for (j in sel.snps) # for each snp get kinship
      K <- K + geno$probs[,,j] %*% score.matrix %*% t(geno$probs[,,j])    
      
    # normalization to diag(K) == 1
    norm.factor <- sqrt(diag(K))
    K <- diag(1/norm.factor) %*% K %*%  diag(1/norm.factor)
    
  }
  
  if (method == 'allele-multif-additive') {
    score.matrix <- diag(length(geno$calls))
    K <- matrix(0, nrow(geno$probs), nrow(geno$probs))
    for (j in sel.snps) # for each snp get kinship
      K <- K + geno$probs[,,j] %*% score.matrix %*% t(geno$probs[,,j])
      
    # normalization to diag(K) == 1
    norm.factor <- sqrt(diag(K))
    K <- diag(1/norm.factor) %*% K %*%  diag(1/norm.factor)
    
  }
  
  if (method == 'allele-multif-cosine') {
    
    K <- matrix(0, nrow(geno$probs), nrow(geno$probs))
    diag(K) <- 1
    
    for (i in seq(nrow(geno$probs)-1))
      for (j in (i+1):nrow(geno$probs))
        K[i,j] <- K[j,i] <- sum(geno$probs[i,,sel.snps] * geno$probs[j,,sel.snps]) / sqrt( sum(geno$probs[i,,sel.snps]^2) * sum(geno$probs[j,,sel.snps]^2) )
  }
  
  return(K)
  
}