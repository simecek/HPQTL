#' Check \code{genotype.probs} object
#'  
#' @param geno a \code{genotype.probs} object
#'
#' @return \code{NULL}
#' 
#' @keywords manip
#'
#' @examples
#' data(fake.f2, package="qtl")
#' fake.f2 <- calc.genoprob(fake.f2)
#' 
#' geno <- extract.geno(fake.f2)
#' HPQTL:::check.genotype.probs

check.genotype.probs <- function(geno, lmm.l1o=FALSE) {
  stopifnot("genotype.probs" %in% class(geno)) # geno is genotype.probs" class
  stopifnot(class(geno$probs) == "array")
  stopifnot( dim(geno$probs)[1] == length(geno$subjects) ) # number of subjets equals
  stopifnot( dim(geno$probs)[2] == length(geno$calls) ) # number of possible calls equals
  stopifnot( dim(geno$probs)[3] == nrow(geno$markers) ) # number of markers equals
  if (lmm.l1o) 
  if (!is.null(geno$chromosomes)) stopifnot(sort(unique(geno$markers$chr)) == sort(geno$chromosomes$chr)) 
}