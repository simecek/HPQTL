#' Add a QTL to Simulated Phenotype
#' 
#' @param cross "\code{cross}" object
#' @param pheno.cols selection of phenotype's column(s)
#' @param qtl.pos position of QTL(s), either number of markers or relative position in genome (between 0 and 1)
#' @param qtl.size size of QTL(s), 
#'
#' @return "\code{cross}" object
#' 
#' @keywords manip
#'
#' @export
#' 
#' @examples
#' cross <- sim.cross.geno(250, nmar=10)
#' cross$pheno <- sim.cross.pheno(0.01, cross)
#' plot(scan1(cross))
#' cross <- add.qtl(cross, 0.75, 1)
#' plot(scan1(cross))

add.qtl <- function(cross, qtl.pos, qtl.size, pheno.cols=1) {
  
  stopifnot(length(qtl.pos) == length(qtl.size))
  geno <- extract.geno(cross)
  
  # if position is relative, find the marker index
  if (all(qtl.pos<1)) qtl.pos <- floor(qtl.pos * ncol(geno)) + 1
  
  for (p in pheno.cols) {
    #sd.geno <- apply(geno[, qtl.pos, drop=FALSE], 2, sd)
    sd.pheno <- sd(cross$pheno[,p])
    cross$pheno[,p] <- cross$pheno[,p] + sd.pheno * colSums(t(geno[, qtl.pos, drop=FALSE] - 2) * sqrt(qtl.size))
  }
  
  cross
}

