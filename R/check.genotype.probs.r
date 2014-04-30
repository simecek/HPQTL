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
  stopifnot(all.equal(length(geno$probs), length(geno$markers), length(geno$calls))) # same length of lists
  stopifnot( length(unique(sapply(geno$probs, dim)[1,])) == 1 ) # number of subjets is equal
  stopifnot( sapply(geno$probs, dim)[2,] == sapply(geno$calls,length) ) # number of possible calls match
  stopifnot( sapply(geno$probs, dim)[3,] == sapply(geno$markers,nrow) ) # number of possible markers match
}