#' Genome Scan With a Single QTL Model
#' 
#' @param geno genotype probabilities ("\code{genotype.probs}" or "\code{cross}" object)
#' @param pheno data frame with phenotypes
#' @param pheno.cols selection of phenotype's column(s)
#' @param covar (additive) covariates
#' @param intcovar interactive covariates (interact with QTL genotype)
#' @param procedure procedure to do the inference, see Details
#' @param G genetic similarity matrix or a list of genetic similarity matrices or \code{NULL}
#' @param Intercept option to pass rotated intercept for LMM model
#' @param ... parameters passed to \code{gensim.matrix}
#' 
#' @details Currently, three procedures are implemented: linear model (LM), linear mixed model (LMM)
#' and linear mixed model - leave one chromosome out (LOCO).
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
#' plot(scan1(fake.f2, procedure="LOCO"), incl.markers=FALSE)

scan1 <- function(geno, pheno, pheno.col=1, rankZ=FALSE, covar=NULL, intcovar=NULL,
                  procedure=c("LM","LMM","LOCO"), G, verbose=FALSE,
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
    G <- gensim.matrix(geno, procedure=procedure) 
  } else {
    if (procedure == "LM" & !is.null(G)) warning("For 'LM' procedure, gen. sim. matrix G is not used")
    if (procedure == "LMM" & !is.matrix(G)) stop("For 'LMM' procedure, G should be matrix")
    if (procedure == "LOCO" & !is.list(G)) stop("For 'LOCO' procedure, G should be a list of matrices")
  }  
  
  # complete cases only
  if (rankZ) y <- HPQTL:::rankZ(pheno[,pheno.col]) else y <- pheno[,pheno.col]
  selected <- which(complete.cases(cbind(y,covar)))
  n.selected <- length(selected) # number of individuals
  if (is.null(covar)) covar <- cbind(Intercept)
  
  # prepare output
  output <- data.frame(chr=as.character(geno$markers$chr), 
                       pos=geno$markers$pos, 
                       lod=matrix(0, nrow=nrow(geno$markers), ncol = 1), 
                       row.names=geno$markers$marker,
                       stringsAsFactors = FALSE)
  class(output) <- c("scanone", "data.frame")
  
  if (verbose) t1 <- Sys.time() # measure time for variance decomposition
  
  if (procedure == "LMM") {
    if (verbose) message("Variance decomposition started.")
    A <- variance.decomposition(y[selected], covar[selected,], G[selected,selected],...)$A
  }
  if (procedure == "LOCO") {
    A <- foreach (c = geno$chromosomes$chr) %do% {
      if (verbose) message(paste("Variance decomposition for chromosome", c))
      variance.decomposition(y[selected], covar[selected,], G[[c]][selected,selected],...)$A
    }
    names(A) <- geno$chromosomes$chr
  }
  
  if (verbose) t2 <- Sys.time() # measure time for variance decomposition
  if (verbose) message(paste("Variance decomposition takes", t2-t1, units(t2-t1)))
  
  ### Rotate probs, y and covars

  
  if (verbose) t1 <- Sys.time() # measure time for rotation
  
  # rotate genotype probabilities
  probs.rot <- foreach (c = geno$chromosomes$chr) %do% {
    if (verbose) message(paste("Rotating chromosome", c))
    
    switch(procedure,
           LM = geno$probs[selected,-1,geno$markers$chr==c, drop=FALSE],
           LMM = { 
             tmp <- geno$probs[selected,-1,geno$markers$chr==c, drop=FALSE]
             dimtmp <- dim(tmp) 
             dim(tmp) <- c(dim(tmp)[1], dim(tmp)[2]*dim(tmp)[3])
             tmp <- A %*% tmp
             dim(tmp) <- dimtmp
             tmp
           },
           LOCO = { 
             tmp <- geno$probs[selected,-1,geno$markers$chr==c, drop=FALSE]
             dimtmp <- dim(tmp) 
             dim(tmp) <- c(dim(tmp)[1], dim(tmp)[2]*dim(tmp)[3])
             tmp <- A[[c]] %*% tmp
             dim(tmp) <- dimtmp
             tmp
           })
  }
  names(probs.rot) <- geno$chromosomes$chr

  if (!is.null(intcovar)) {
    # rotate genotype probabilities
    inter.rot <- foreach (c = geno$chromosomes$chr) %do% {
      if (verbose) message(paste("Interactive covariate for chromosome", c))
      
      switch(procedure,
             LM = geno$probs[selected,-1,geno$markers$chr==c, drop=FALSE] * intcovar[selected],
             LMM = { 
               tmp <- geno$probs[selected,-1,geno$markers$chr==c, drop=FALSE]
               dimtmp <- dim(tmp) 
               dim(tmp) <- c(dim(tmp)[1], dim(tmp)[2]*dim(tmp)[3])
               tmp <- tmp * intcovar[selected]
               tmp <- A %*% tmp
               dim(tmp) <- dimtmp
               tmp
             },
             LOCO = { 
               tmp <- geno$probs[selected,-1,geno$markers$chr==c, drop=FALSE]
               dimtmp <- dim(tmp) 
               dim(tmp) <- c(dim(tmp)[1], dim(tmp)[2]*dim(tmp)[3])
               tmp <- tmp * intcovar[selected]
               tmp <- A[[c]] %*% tmp
               dim(tmp) <- dimtmp
               tmp
             })
    }
    names(inter.rot) <- geno$chromosomes$chr
  }
  
  # rotate y
  y.rot <- foreach (c = geno$chromosomes$chr) %do% {
    switch(procedure,
           LM = cbind(y[selected]),
           LMM = A %*% cbind(y[selected]),
           LOCO = A[[c]] %*% cbind(y[selected]))
  }
  names(y.rot) <- geno$chromosomes$chr
  
  covar.rot <- foreach (c = geno$chromosomes$chr) %do% {
    switch(procedure,
           LM = covar[selected,],
           LMM = A %*% covar[selected,],
           LOCO = A[[c]] %*% covar[selected,])
  }
  names(covar.rot) <- geno$chromosomes$chr
  
  if (verbose) t2 <- Sys.time() # measure time for variance decomposition
  if (verbose) message(paste("Rotation takes", t2-t1, units(t2-t1)))
 
  if (verbose) t1 <- Sys.time() # measure time for LOD calculation
  
  # null models
  rss0 <- foreach (c = geno$chromosomes$chr, .combine="c") %do% {
    tmp <- sum(lsfit(y=y.rot[[c]], x=covar.rot[[c]], intercept=FALSE)$residuals^2)
    rep(tmp, dim(probs.rot[[c]])[3])
  }
  
  # model with additive QTL
  rss1 <- foreach (c = geno$chromosomes$chr, .combine="c") %do% {
    if (verbose) message(paste("Processing chromosome", c))
    foreach (i = 1:dim(probs.rot[[c]])[3], .combine="c") %do% {     
      sum(lsfit(y=y.rot[[c]], x=cbind(covar.rot[[c]],probs.rot[[c]][,,i]), intercept=FALSE)$residuals^2)
    }
  }  
  
  output$lod <-  n.selected/2 * (log10(rss0) - log10(rss1))
  
  # model with intcovar
  if (!is.null(intcovar)) {
    w <- getOption("warn")
    options(warn = -1)
    rss2 <- foreach (c = geno$chromosomes$chr, .combine="c") %do% {
      if (verbose) message(paste("Processing chromosome", c))
      foreach (i = 1:dim(probs.rot[[c]])[3], .combine="c") %do% {     
        sum(lsfit(y=y.rot[[c]], x=cbind(covar.rot[[c]], probs.rot[[c]][,,i], inter.rot[[c]][,,i]), intercept=FALSE)$residuals^2)
      }
    }
    options(warn = w)
    
    output$lod2 <-  n.selected/2 * (log10(rss1) - log10(rss2))
    output$lod3 <-  n.selected/2 * (log10(rss0) - log10(rss2))
  }
  
  if (verbose) t2 <- Sys.time() # measure time for variance decomposition
  if (verbose) message(paste("LOD calculation takes", t2-t1, units(t2-t1)))

  if (verbose) message("Finished")
  
  return(output)
}
