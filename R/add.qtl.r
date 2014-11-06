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
#' geno <- sim.cross.geno(250, nmar=10)
#' pheno <- sim.cross.pheno(0.01, geno)
#' plot(scan1(geno, pheno))
#' pheno2 <- add.qtl(geno, pheno, 0.75, 1)
#' plot(scan1(geno, pheno2))

add.qtl <- function(geno, pheno, qtl.pos, qtl.size, pheno.cols=1) {
  
  stopifnot(length(qtl.pos) == length(qtl.size))
  
  # if position is relative, find the marker index
  if (all(qtl.pos<1)) qtl.pos <- floor(qtl.pos * length(geno$markers)) + 1
  
  for (p in pheno.cols) {
    sd.geno <- as.vector(apply(geno$probs[,1,qtl.pos,drop=FALSE],c(2,3),sd))
    sd.pheno <- sd(pheno[,p])
    pheno[,p] <- pheno[,p] + sd.pheno * apply(aperm(geno$probs[, 1,qtl.pos, drop=FALSE], c(3,2,1)) * sqrt(qtl.size)/sd.geno, c(2,3), sum)
  }
  
  pheno
}

