#' Genome Scan With a Single QTL Model
#' 
#' @param geno genotype probabilities ("\code{genotype.probs}" or "\code{cross}" object)
#' @param pheno data frame with phenotypes
#' @param pheno.cols selection of phenotype's column(s)
#' @param covar (additive) covariates
#' @param procedure procedure to do the inference, see Details
#' @param G genetic similarity matrix or a list of genetic similarity matrices or \code{NULL}
#' @param subjects subseting of subjects
#' @param markers subseting of markers
#' @param Intercept option to pass rotated intercept for LMM model
#' @param ... parameters passed to \code{gensim.matrix}
#' 
#' @details Currently, three procedures are implemented: linear model (LM), linear mixed model (LMM)
#' and linear mixed model - leave one chromosome out (LMM-L1O).
#' 
#' @return \code{scanone} object
#' 
#' @keywords manip
#'
#' @export
#' 
#' @examples
#' data(fake.f2, package="qtl")
#' fake.f2 <- calc.genoprob(fake.f2)
#' 
#' plot(scan1(fake.f2, procedure="LM"), incl.markers=FALSE)
#' plot(scan1(fake.f2, procedure="LMM"), incl.markers=FALSE)
#' plot(scan1(fake.f2, procedure="LMM-L1O"), incl.markers=FALSE)

scan1 <- function(geno, pheno, pheno.cols=1, covar=NULL, procedure=c("LM","LMM","LMM-L1O"), G, 
                  subjects=seq(geno$subjects), markers=seq(NROW(geno$markers)), 
                  Intercept=rep(1,length(geno$subjects)), ...) {
 
  procedure <- match.arg(procedure)

  # covar should be matrix, not data.frame
  if (!is.null(covar) & class(covar)!="matrix") {
    warning(paste("Class of 'covar' should be matrix, not", class(covar)))
    covar = as.matrix(covar)
  }
  
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
  
  # prepare output
  output <- data.frame(chr=geno$markers$chr, 
                       pos=geno$markers$pos, 
                       lod=matrix(0, nrow=nrow(geno$markers), ncol = length(pheno.cols)), 
                       row.names=geno$markers$marker)
  class(output) <- c("scanone", "data.frame")
  
  for (i in pheno.cols) {
    
    y <- pheno[,i]
    index.i = which(pheno.cols==i)[1]
    selected <- intersect(subjects, which(complete.cases(cbind(y,covar))))  
    n.selected <- length(selected) # number of individuals

    
    # linear model
    if (procedure == "LM") A <- NULL
    # linear mixed model
    if (procedure == "LMM")
      A <- variance.decomposition(y[selected], covar[selected,], G[selected,selected],...)$A
    # linear mixed model, leave one out
    if (procedure == "LMM-L1O") {
      A <- list()
      for (c in geno$chromosomes$chr) {
        A[[c]] <- variance.decomposition(y[selected], covar[selected,], G[[c]][selected,selected],...)$A
      }
    } 
    
    # LOD caclulations are performed per chromosome
    for (c in geno$chromosomes$chr) {
      
      # found appropriate A.chr matrix
      if (procedure == "LMM-L1O") A.chr <- A[[c]] else A.chr <- A
      
      # unless LM procedure, rotate y and covar
      if (procedure != "LM") {
        y.rotated <- A.chr %*% y[selected] 
        covar.rotated <- A.chr %*% cbind(Intercept,covar)[selected,]
      } else {
        y.rotated <- y[selected] 
        covar.rotated <- cbind(Intercept,covar)[selected,]
      }
      
      # null model
      rss0 <- sum(lsfit(y=y.rotated, x=covar.rotated, intercept=FALSE)$residuals^2)
      
      for (j in intersect(which(geno$markers$chr == c),markers)){
        if (procedure != "LM") 
          x.rotated <- A.chr %*% cbind(geno$probs[selected,,j],covar[selected,]) 
        else 
          x.rotated <- cbind(geno$probs[selected,,j], covar[selected,])
        
        rss1 <- sum(lsfit(y=y.rotated, x=x.rotated, intercept=FALSE)$residuals^2)
        output[j,index.i+2] <- n.selected/2 * (log10(rss0) - log10(rss1))
      } 
    }
  }
    
  return(output[markers,]) 
  
}