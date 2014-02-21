#' Extract Genotype Matrix
#' 
#' @param cross "\code{cross}" object
#'
#' @return matrix (\code{nind} rows, \code{totalmar} columns)
#' 
#' @keywords manip
#'
#' @export
#' 
#' @examples
#' cross <- sim.cross.geno(250, nmar=10)
#' extract.geno(cross)

extract.geno <- function(cross) {
  tmp <- lapply(cross$geno, function(x) x[[1]])
  do.call("cbind", tmp)
}