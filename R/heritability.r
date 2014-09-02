#' Estimate Heritability
#' 
#' @param geno genotype probabilities ("\code{genotype.probs}" or "\code{cross}" object)
#' @param pheno data frame with phenotypes
#' @param pheno.cols selection of phenotype's column(s)
#' @param covar (additive) covariates
#' @param G genetic similarity matrix or a list of genetic similarity matrices or \code{NULL}
#' @param se if se=TRUE standard error is estimates
#' @param ... parameters passed to \code{gensim.matrix}
#' @return numeric
#' 
#' @keywords manip
#'
#' @export

heritability <- function(geno, pheno, pheno.cols=1, covar=NULL, se=FALSE, G, ...) {
  
  # if geno is not 'genotype.probs', export genotype
  if (!('genotype.probs' %in% class(geno))) {
    # if phenotype is missing, try to extract it
    if (missing(pheno) & "pheno" %in% names(geno))
      pheno <- geno$pheno
    geno <- extract.geno(geno)
  } 
  
  # if G is missing then it should be estimated from genotype  
  if (missing(G)) {
    G <- gensim.matrix(geno, ...)
  }
    
  output <- c()
  if (se) output.se <- c()
  
  for (p in pheno.cols) {
  
    # fit the mixed model
    y <- pheno[,p]
    rg.fit <- regress(y~covar, ~G)
  
    # estimate heritability
    h2 <- as.numeric(rg.fit$sigma[1] / sum(rg.fit$sigma))
    output <- c(output, h2)
    if (se) {
      v <- c(rg.fit$sigma[2], -rg.fit$sigma[1]) / sum(rg.fit$sigma^2)
      h2.se <- rbind(v) %*% rg.fit$sigma.cov %*% cbind(v)
      output.se <- c(output.se, h2.se)
    }  
  }
  
  if (se) attr(output, "se") <- sqrt(output.se)
  output
}