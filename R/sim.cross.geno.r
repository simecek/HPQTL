#' Generate a Genotype of Random Cross
#' 
#' @param nind number of animals
#' @param ctype cross type ("f2" / "bc")
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

sim.cross.geno <- function(nind = 250, ctype = c("f2", "bc"), nchr=20, chrlen=100, nmar=10) {
  ctype <- match.arg(ctype)
  map <- sim.map(len=rep(chrlen, nchr), n.mar=nmar)
  cross <- sim.cross(map, type=ctype, n.ind=nind, model = NULL)
  cross <- calc.genoprob(cross)
  
  # supress warning because of chrX issue of extract.geno
  op <- options()
  options(warn=-1)
  output <- extract.geno(cross)
  options(op)
  output
}