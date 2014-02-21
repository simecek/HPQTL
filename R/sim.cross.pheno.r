#' Generate Random Phenotype with Given Heritability
#' 
#' @param h2 heritability
#' @param cross \code{cross} object
#' @param background method to generate genetic background
#' @param ... parameters passed to \code{genrel.matrix}
#'
#' @return "\code{cross}" object
#' 
#' @keywords manip
#'
#' @export
#' 
#' @examples
#' cross <- sim.cross.geno(250, nmar=10)
#' cross$pheno <- sim.cross.pheno(0.5, cross)

sim.cross.pheno <- function(h2, cross, background = "GRM", ...) {
  # extract genotype matrix
  geno <- extract.geno(cross)
  
  if (background == "GRM") {
    G <- genrel.matrix(geno, ...)
    genetic <- mvrnorm(n = 1, rep(0,nind(cross)), G)
    noise <- rnorm(nind(cross))
  }  
  
  if (background == "all-snps") {
    genetic <- colMeans(t(geno - 2) * rnorm(ncol(geno)))
    noise <- rnorm(nind(cross), sd = sd(genetic))
  }  

  output <- NULL
  for (i in 1:length(h2))
    output <- cbind(output, genetic + sqrt(1/h2[i]-1)*noise)
  colnames(output) <- paste0("h", h2)
  return(output)
  
}