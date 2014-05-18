#' Significance Threshold of Genome Scan
#'  
#' @param cross "\code{cross}" object
#' @param pheno.cols selection of phenotype's column(s)
#' @param procedure procedure to do the inference, see Details
#' @param n.perm number of permutations
#' @param alpha level of significance 
#' @param ... parameters passed to \code{genrel.matrix}
#'
#' @details Currently, three procedures fully implemented: \code{procedure = "qtl"} calls \code{scanone} 
#' function from \code{qtl} package, \code{procedure = "QTLRel"} calls \code{scanOne} function from \code{QTLRel} package,
#' \code{procedure = "QTLRel-pre-chr"} calls \code{scanOne} for each chromosome separately and genetic relationship matrix
#' is estimated with this chromosome excluded.
#'
#' @return numeric
#' 
#' @keywords manip
#'
#' @export
#' 
#' @examples
#' cross <- sim.cross.geno(250, nmar=10)
#' cross$pheno <- sim.cross.pheno(0.5, cross)
#' scan1.threshold(cross)
#' 
#' # slow - do not run
#' # scan1.threshold(cross, n.perm=1000, alpha=0.01)
#' # scan1.threshold(cross, procedure = "scanOne")
#' # scan1.threshold(cross, procedure = "scanOne-per-chr")

scan1.threshold <- function(pheno, geno, covar=NULL, procedure=c("LM","LMM","LMM-L1O"), G=NULL, n.perm=100, alpha=0.05, ...) {
  
  procedure <- match.arg(procedure)
  nonNA <- which(!is.na(pheno))
  
  if (procedure == "LM") {
    
    maxlods <- replicate(n.perm, max(scan1(pheno, geno, covar, procedure="LM", NULL, shuffle=TRUE)$lod))
    trhold <- quantile(maxlods, 1 - alpha)
  }
    
  if (procedure == "LMM") {
    
    if (is.null(G)) G <- gensim.matrix(geno, ...)
    vd <- variance.decomposition(pheno[nonNA], covar[nonNA,], G[nonNA,nonNA],...)
    
    Intercept <- rep(1, length(nonNA))
    Yperm <- mvnpermute(pheno[nonNA], cbind(Intercept, covar[nonNA,]), vd$V, n.perm)    
    Yperm.rotated <- vd$A %*% Yperm
    if (!is.null(covar)) covar.rotated <- vd$A %*% covar[nonNA,] else covar.rotated <- NULL
    Intercept.rotated = vd$A %*% rep(1, length(nonNA))
    maxlods <- rep(0, n.perm)
    
    for (i in seq(n.perm))
      maxlods[i] <- max(scan1(Yperm.rotated[,i], geno, covar.rotated, procedure="LM", NULL, Intercept = Intercept.rotated)$lod)
     
    trhold <- quantile(maxlods, 1 - alpha)
  }
  
  if (procedure == "LMM-L1O") {
      
    if (is.null(geno$chromosomes) | is.null(geno$markers$chr) | is.null(geno$markers$pos))
      stop("For LMM-L1O procedure, markers MUST be mapped to chromosomes.")
    if (is.null(G)) {
      warning("Genetic similarity matrix G should be specified, use 'gensim.matrix' function.")
      G <- list()
      for (c in geno$chromosomes$chr)
        G[[c]] <- gensim.matrix(geno, skip.chr=c, ...)
    }
    
    maxlods <- matrix(nrow = NROW(geno$chromosomes), ncol = n.perm)

    for (c in geno$chromosomes$chr) {
        vd <- variance.decomposition(pheno[nonNA], covar[nonNA,], G[[c]][nonNA,nonNA],...)
        
        Y.perm <- mvnpermute(pheno[nonNA], cbind(Intercept, covar[nonNA,]), vd$V, n.perm)
        Yperm.rotated <- vd$A %*% Yperm
        if (!is.null(covar)) covar.rotated <- vd$A %*% covar[nonNA,] else covar.rotated <- NULL
        Intercept.rotated = vd$A %*% rep(1, length(nonNA))
        
        for (i in seq(n.perm))
          maxlods[which(chrs==c),i] <- max(scan1(Yperm.rotated[,i], geno[,,geno$markers$chr==c], covar.rotated, procedure="LM", NULL, Intercept = Intercept.rotated)$lod)
    }

    trhold <- quantile(apply(maxlods, 2, max), 1-alpha)
  }  
  
   
  return(as.numeric(trhold))
 
}