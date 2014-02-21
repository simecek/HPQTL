#' Estimate Heritability
#' 
#' @param cross "\code{cross}" object
#' @param pheno.cols selection of phenotype's column(s)
#' @param ... parameters passed to \code{genrel.matrix}
#'
#' @return numeric
#' 
#' @keywords manip
#'
#' @export
#' 
#' @examples
#' cross <- sim.cross.geno(250, nmar=10)
#' cross$pheno <- sim.cross.pheno(0.5, cross)
#' heritability(cross)

heritability <- function(cross, pheno.cols=1, ...) {
  
  # get genetic relationship matrix
  G <- genrel.matrix(cross, ...)
  output <- c()
  
  for (p in pheno.cols) {
  
    # fit the mixed model
    Y <- cross$pheno[,p]
    rg.fit <- regress(Y~1, ~G)
  
    # estimate heritability
    h2 <- as.numeric(rg.fit$sigma[1] / sum(rg.fit$sigma))
    output <- c(output, h2)
  }
  
  output
}