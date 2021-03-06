#' Generate Random Phenotype with Given Heritability
#' 
#' @param h2 heritability
#' @param geno "\code{genotype.probs}" object
#' @param background method to generate genetic background, see Details
#' @param ... parameters passed to \code{genrel.matrix}
#'
#' @details Currently two polygenic backgrounds are supported: 
#' If \code{background = "GRM"} then the genetic effect is generated from a multivariate 
#' normal distribution with genetic similarity matrix as variance matrix. If 
#' \code{background = "all-snps"} then genetic effect is a sum of small gaussian distributed
#' effect at every SNP location.
#' 
#' @return matrix with phenotype data
#' 
#' @keywords manip
#'
#' @export
#' 
#' @examples
#' cross <- sim.cross.geno(250, nmar=10)
#' cross$pheno <- sim.cross.pheno(0.5, cross)
#' cross$pheno <- sim.cross.pheno(0.5, cross, method = "kinship")
#' cross$pheno <- sim.cross.pheno(0.5, cross, background = "all-snps")

sim.cross.pheno <- function(h2, geno, background = "GRM", ...) {
  
  n <- length(geno$subjects)
  
  if (background == "GRM") {
    G <- gensim.matrix(geno, ...)
    genetic <- mvrnorm(n = 1, rep(0,n), G)
    noise <- rnorm(n)
  }  
  
  output <- NULL
  for (i in 1:length(h2))
    output <- cbind(output, genetic + sqrt(1/h2[i]-1)*noise)
  colnames(output) <- paste0("h", h2)
  return(output)
  
}