#' Generate a Genotype of Random Cross
#' 
#' @param nind number of animals
#' @param ctype cross type ("f2" / "bc" / "risib")
#' @param nchr number of chromosomes
#' @param chrlen length of chromosomes (in cM)
#' @param nmar number of markers per chromosome
#'
#' @return "\code{cross}" object
#' 
#' @keywords manip
#'
#' @export
#' 
#' @examples
#' sim.cross.geno(250, nmar=10)

sim.cross.geno <- function(nind = 250, ctype = "f2", nchr=20, chrlen=100, nmar=250) {
  map <- sim.map(len=rep(chrlen, nchr), n.mar=nmar)
  sim.cross(map, type=ctype, n.ind=nind, model = NULL)
}