#' Add a QTL to Simulated Phenotype
#' 
#' @param cross "\code{cross}" object
#' @param pheno.cols selection of phenotype's column(s)
#' @param qtl.pos position of QTL(s), qither number of markets or relative position in genome (between 0 and 1)
#' @param qtl.size size of QTL(s), 
#'
#' @return numeric
#' 
#' @keywords manip
#'
#' @export
#' 
#' @examples
#' cross <- sim.cross.geno(250, nmar=10)
#' cross <- calc.genoprob(cross)
#' plot(scanone(cross))
#' cross.with.qtls <- add.qtl(cross, c(0.15,0.35), c(1,1))
#' plot(scanone(cross.with.qtls))

add.qtl <- function(cross, qtl.pos, qtl.size, pheno.cols=1) {
  
  stopifnot(length(qtl.pos) == length(qtl.size))
  
  geno <- extract.geno(cross)
  
  if (all(qtl.pos<1)) qtl.pos <- floor(qtl.pos * ncol(geno)) + 1
  pheno.sd <- apply(cross$pheno, 2, sd)
  
  for (p in pheno.cols) {
    sd.geno <- apply(geno[, qtl.pos, drop=FALSE], 2, sd)
    sd.pheno <- sd(cross$pheno[,p])
    cross$pheno[,p] <- cross$pheno[,p] + sd.pheno * colSums(t(geno[, qtl.pos, drop=FALSE] - 2) * sqrt(qtl.size) / sd.geno)
  }
  
  cross
}

