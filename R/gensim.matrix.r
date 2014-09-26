#' Genetic Similarity Matrix
#' 
#'   
#' @param geno genotype probabilities ("\code{genotype.probs}" or "\code{cross}" object)
#' @param method how to calculate similarity matrix
#' @param procedure procedure of scan1 the G is calculated for
#' @param subjects subseting of subjects
#' @param markers subseting of markers
#' @param gensim.normalization method to be used for genetic similarity matrix normalization   
#'
#' @return square matrix (\code{nind} rows, \code{nind} columns)
#' 
#' @keywords manip
#'
#' @export
#' 
#' @examples
#' data(fake.f2, package="qtl")
#' fake.f2 <- calc.genoprob(fake.f2)
#' 
#' geno <- extract.geno(fake.f2)
#' G <- gensim.matrix(geno)
#' Glist <- gensim.matrix(geno, procedure="LOCO")

gensim.matrix <- function(geno, method=c('default', 'allele-2f-additive', 'allele-multif-additive', 'allele-multif-cosine'),
                          procedure = c("LMM", "LOCO", "LM"), subjects=seq(geno$subjects), markers = seq(NROW(geno$markers)), 
                          gensim.normalization="sample-variance", ...) {
  
  # match arguments
  method = match.arg(method) 
  
  procedure = match.arg(procedure)
  
  # if not genotype.probs then extract genotype probs
  if ("cross" %in% class(geno)) geno <- extract.geno(geno)
  check.genotype.probs(geno)

  # for default 'method' guess 
  if (method=="default")
    if (length(geno$calls)<=3) 
      method <- 'allele-2f-additive'
    else 
      method <- 'allele-multif-cosine'
    
  # for linear model no genetic similarity matrix is needed
  if (procedure == "LM") return(NULL)
  
  # for LOCO call it recursively for each cheomosome
  if (procedure == "LOCO") {
    Glist <- list()
    for (c in geno$chromosomes$chr)
      Glist[[c]] <- gensim.matrix(geno, method=method, procedure="LMM", subjects=subjects, 
                                  markers = intersect(markers,which(geno$markers$chr!=c)), 
                                  gensim.normalization=gensim.normalization)
    return(Glist)
  }
  
  if (method == 'allele-2f-additive') {
    if (length(geno$calls)<2 | length(geno$calls)>3) stop("Method 'allele-2f-additive' expects 2 founders.")
    score.matrix <- switch(length(geno$calls),
                           NULL,
                           matrix(c(1,1/2,1/2,1),2,2),                     
                           matrix(c(1,1/2,0,1/2,1,1/2,0,1/2,1),3,3)) 
    K <- matrix(0, length(subjects), length(subjects))
    
    for (j in markers) # for each snp get kinship
      K <- K + geno$probs[subjects,,j] %*% score.matrix %*% t(geno$probs[subjects,,j])    
    
  }
  
  if (method == 'allele-multif-additive') {
    score.matrix <- diag(length(geno$calls))
    K <- matrix(0, length(subjects), length(subjects))
    for (j in markers) # for each snp get kinship
      K <- K + geno$probs[subjects,,j] %*% score.matrix %*% t(geno$probs[subjects,,j])
    
  }
  
  if (method == 'allele-multif-cosine') {
    
    mprobs = geno$probs[subjects,,markers]
    dp = dim(mprobs)
    dim(mprobs) = c(dp[1], dp[2]*dp[3])
    
    K <- mprobs %*% t(mprobs)
  }
  
  # normalize the matrix
  K <- normalize.matrix(K, method=gensim.normalization)
  
  return(K)
  
}