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

gensim.matrix <- function(geno, method=c('default', 'allele-2f-additive', 'allele-multif-additive', 'allele-multif-cosine'),
                          procedure = c("LMM", "LMM-L1O", "LM"), subjects=seq(geno$subjects), markers = seq(NROW(geno$markers))) {
  
  # match arguments
  method = match.arg(method)
  procedure = match.arg(procedure)
  
  # if not genotype.probs then extract genotype probs
  if ("cross" %in% class(geno)) geno <- extract.geno(geno)
  check.genotype.probs(geno)

  # for linear model no genetic similarity matrix is needed
  if (procedure == "LM") G <- NULL
  
  # for LMM-L1O call it recursively for each cheomosome
  if (procedure == "LMM-L1O") {
    G <- list()
    for (c in geno$chromosomes$chr)
      G[[c]] <- gensim.matrix(geno, method=method, procedure="LMM", subjects=subjects, markers = which(geno$markers$chr))
  }
  
  if (method == 'allele-2f-additive') {
    if (length(geno$calls)<2 | length(geno$calls)>3) stop("Method 'allele-2f-additive' expects 2 founders.")
    score.matrix <- switch(length(geno$calls),
                           NULL,
                           matrix(c(1,1/2,1/2,1),2,2),                     
                           matrix(c(1,1/2,0,1/2,1,1/2,0,1/2,1),3,3)) 
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