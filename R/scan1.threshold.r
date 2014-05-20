#' Significance Threshold of Genome Scan
#'  
#' @param geno genotype probabilities ("\code{genotype.probs}" or "\code{cross}" object)
#' @param pheno data frame with phenotypes
#' @param pheno.cols selection of phenotype's column(s)
#' @param covar (additive) covariates
#' @param procedure procedure to do the inference, see Details
#' @param G genetic similarity matrix or a list of genetic similarity matrices or \code{NULL}
#' @param subjects subseting of subjects
#' @param markers subseting of markers
#' @param n.perm number of permutations
#' @param alpha level of significance 
#' @param ... parameters passed to \code{genrel.matrix}
#'
#' @details Currently, three procedures are implemented: linear model (LM), linear mixed model (LMM)
#' and linear mixed model - leave one chromosome out (LMM-L1O).
#'
#' @return numeric
#' 
#' @keywords manip
#'
#' @export
#' 
#' @examples
#' data(fake.f2, package="qtl")
#' fake.f2 <- calc.genoprob(fake.f2)
#' 
#' # warning, n.perm=10 is too low for practical purposes (but fast)
#' scan1.threshold(fake.f2, procedure="LM", n.perm=10)
#' scan1.threshold(fake.f2, procedure="LMM", n.perm=10)
#' scan1.threshold(fake.f2, procedure="LMM-L1O", n.perm=10)

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