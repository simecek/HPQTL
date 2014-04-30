#' Extract Genotype Probabilities
#'
#' @details Extract genotype prob. array from "\code{qtl::cross}" object. 
#' Stored as a list  
#'   
#' @param cross "\code{cross}" object
#'
#' @return \code{genotype.probs} object is a list of three components 
#' \itemize{ 
#'   \item \code{probs} 3-dimensional array of genotype probabilities (subjects x calls x markers) or a list of such arrays
#'   \item \code{calls} possible genotype calls 
#'   \item \code{markers}
#' } 
#' 
#' @keywords manip
#'
#' @export
#' 
#' @examples
#' cross <- sim.cross.geno(250, nmar=10)
#' extract.geno(cross)

extract.geno <- function(cross) {

  #currently only cross (qtl package) implemented
  stopifnot("cross" %in% class(cross))
  
  if ("cross" %in% class(cross)) {
  
    # check that calc.genoprob has been called
    if (!("prob" %in% names(cross$geno[[1]]))) {
      warning("First running calc.genoprob.")
      cross <- calc.genoprob(cross)
    }
    
    # genotype is a 3dim array of probs (subjects x calls x markers), vector of calls and vector of markers 
    # or lists of such triples
    geno <- list(probs=list(), calls=list(), markers=list())  
    for (i in 1:length(cross$geno)) {
      geno$probs[[i]] <- aperm(cross$geno[[i]]$prob, c(1,3,2))
      geno$calls[[i]] <- dimnames(cross$geno[[i]]$prob)[[3]]
      geno$markers[[i]] <- data.frame(marker = names(cross$geno[[i]]$map),
                                      chr = names(cross$geno)[i], 
                                      pos = as.vector(cross$geno[[i]]$map))
    }
    
  }
  
  class(geno) <- c("genotype.probs", "list")
  check.genotype.probs(geno)
  geno
}