#' Check \code{genotype.probs} object
#'  
#' @param geno a \code{genotype.probs} object
#'
#' @return \code{NULL}
#' 
#' @keywords manip
#'
#' @examples
#' cross <- sim.cross.geno(250, nmar=10)
#' geno <- extract.geno(cross)
#' check.genotype.probs(geno)

check.genotype.probs <- function(geno) {
  stopifnot("genotype.probs" %in% class(geno)) # geno is genotype.probs" class
  stopifnot(class(geno$probs) == "array")
  stopifnot( dim(geno$probs)[1] == length(geno$subjects) ) # number of subjets equals
  stopifnot( dim(geno$probs)[2] == length(geno$calls) ) # number of possible calls equals
  stopifnot( dim(geno$probs)[3] == nrow(geno$markers) ) # number of markers equals
  if (!is.null(geno$chromosomes)) stopifnot(sort(unique(geno$markers$chr)) == sort(geno$chromosomes$chr)) 
}