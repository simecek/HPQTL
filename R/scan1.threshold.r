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

scan1.threshold <- function(geno, pheno, pheno.cols=1, covar=NULL, procedure=c("LM","LMM","LMM-L1O"), G, 
                            n.perm=100, alpha=0.05, 
                            subjects=seq(geno$subjects), markers=seq(NROW(geno$markers)), ...) {
  
  procedure <- match.arg(procedure)
  trhold <- c()
  
  # if geno is not 'genotype.probs', export genotype
  if (!('genotype.probs' %in% class(geno))) {
    # if phenotype is missing, try to extract it
    if (missing(pheno) & "pheno" %in% names(geno))
      pheno <- geno$pheno
    geno <- extract.geno(geno)
  }
  
  # if G is missing then it should be estimated from genotype  
  if (missing(G)) {
    G <- gensim.matrix(geno, procedure=procedure, subjects=subjects, markers=markers, ...) 
  } else {
    if (procedure == "LM" & !is.null(G)) warning("For 'LM' procedure, gen. sim. matrix G is not used")
    if (procedure == "LMM" & !is.matrix(G)) stop("For 'LMM' procedure, G should be matrix")
    if (procedure == "LMM-L1O" & !is.list(G)) stop("For 'LMM-L1O' procedure, G should be a list of matrices")
  }
  
  for (i in pheno.cols) {
    
    y <- pheno[,i]
    selected <- intersect(subjects, which(complete.cases(cbind(y,covar))))
    
    if (procedure == "LM") {
      maxlods <- rep(0, n.perm)
      for (j in seq(n.perm)) {
        pheno[selected, i] <- sample(y[selected])
        maxlods[j] <- max(scan1(geno, pheno, pheno.col=i, covar, procedure="LM", 
                                G=NULL, subjects=selected, markers=markers)$lod)
      }  
      trhold <- c(trhold, quantile(maxlods, 1 - alpha))
    }
      
    if (procedure == "LMM") {
      
      vd <- variance.decomposition(y[selected], covar[selected,], G[selected,selected],...)
      
      Intercept = rep(1, length(selected))
      Yperm <- mvnpermute(y[selected], cbind(Intercept, covar[selected,]), vd$V, n.perm)    
      Yperm.rotated <- vd$A %*% Yperm
      if (!is.null(covar)) covar.rotated <- vd$A %*% covar[selected,] else covar.rotated <- NULL
      Intercept.rotated = vd$A %*% Intercept
      for (j in markers)
        geno$probs[selected,,j] <- vd$A %*% geno$probs[selected,,j]
      
      maxlods <- rep(0, n.perm)      
      for (j in seq(n.perm)) {
        pheno[selected,i] = Yperm.rotated[,j]
        maxlods[j] <- max(scan1(geno, pheno, pheno.col=i, covar=covar.rotated, procedure="LM", 
                                G=NULL, Intercept = Intercept.rotated, markers = markers, subjects=selected)$lod)
      }
      
      trhold <- c(trhold, quantile(maxlods, 1 - alpha))
    }
    
    if (procedure == "LMM-L1O") {
        
      maxlods <- matrix(nrow = NROW(geno$chromosomes), ncol = n.perm)
  
      for (c in geno$chromosomes$chr) {
          vd <- variance.decomposition(y[selected], covar[selected,], G[[c]][selected,selected], ...)
          
          Intercept = rep(1, length(selected))
          Yperm <- mvnpermute(y[selected], cbind(Intercept, covar[selected,]), vd$V, n.perm)
          Yperm.rotated <- vd$A %*% Yperm
          if (!is.null(covar)) covar.rotated <- vd$A %*% covar[selected,] else covar.rotated <- NULL
          Intercept.rotated = vd$A %*% Intercept
          cmarkers = intersect(markers, which(geno$markers$chr==c))
          for (j in cmarkers)
            geno$probs[selected,,j] <- vd$A %*% geno$probs[selected,,j]
          
          for (j in seq(n.perm)) {
            pheno[selected,i] = Yperm.rotated[,j]
            maxlods[which(geno$chromosome$chr==c),j] <- max(scan1(geno, pheno, pheno.col=i,  covar=covar.rotated, 
                                                                  procedure="LM", G=NULL, Intercept = Intercept.rotated, 
                                                                  subjects=selected, markers=cmarkers)$lod)
          }  
      }
  
      trhold <- c(trhold,quantile(apply(maxlods, 2, max), 1-alpha))
    }  
  }  
   
  return(as.numeric(trhold))
 
}