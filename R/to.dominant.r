#' Convert 3 calls to 2 (dominance) representation
#'   
#' @param geno genotype probabilities ("\code{genotype.probs}" or "\code{cross}" object) 
#'
#' @return \code{genotype.probs} object with 2 possible calls
#' 
#' @keywords manip
#'
#' 
#' @examples
#' data(fake.f2, package="qtl")
#' fake.f2 <- calc.genoprob(fake.f2)
#' 
#' geno <- extract.geno(fake.f2)
#' geno2 <- to.dominant(geno)


to.additive <- function(geno) {
    
  # if not genotype.probs then extract genotype probs
  if (!("genotype.probs" %in% class(geno))) geno <- extract.geno(geno)
  check.genotype.probs(geno)

  if (length(geno$calls)!=3) stop("Only F2 with 3 calls currently implemented.")
  
  output <- geno
  output$calls <- c("A", "B")
  output$probs <- apply(geno$probs, c(1,3), function(x) c(x[1], x[2]+x[3]))
  output$probs <- aperm(output$probs, c(2,1,3))
  
  output
}