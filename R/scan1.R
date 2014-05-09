#' Genome Scan With a Single QTL Model
#' 
#' @param cross "\code{cross}" object
#' @param pheno.cols selection of phenotype's column(s)
#' @param covar (additive) covariates
#' @param procedure procedure to do the inference, see Details
#' @param V genetic similarity matrix or a list of genetic similarity matrices or \code{NULL}
#' @param ... parameters passed to \code{genrel.matrix}
#' 
#' @details Currently, three procedures fully implemented: \code{procedure = "qtl"} calls \code{scanone} 
#' function from \code{qtl} package, \code{procedure = "QTLRel"} calls \code{scanOne} function from \code{QTLRel} package,
#' \code{procedure = "QTLRel-pre-chr"} calls \code{scanOne} for each chromosome separately and genetic relationship matrix
#' is estimated with this chromosome excluded.
#' 
#' @return \code{scanone} object
#' 
#' @keywords manip
#'
#' @export
#' 
#' @examples
#' cross <- sim.cross.geno(250, nmar=10)
#' cross$pheno <- sim.cross.pheno(c(0.01, 0.5, 0.75), cross)
#' plot(scan1(cross, pheno.cols=1:3), lodcol=1:3)
#' plot(scan1(cross, pheno.cols=1:3, procedure="scanOne"), lodcol=1:3)
#' 
#' # slow - do not run
#' # plot(scan1(cross, pheno.cols=1:3, procedure="scanOne-per-chr"), lodcol=1:3)

scan1 <- function(pheno, geno, covar=NULL, procedure=c("LM","LMM","LMM-L1O"), G=NULL, ...) {

  W.inv<- function(W, symmetric=TRUE,inverse=TRUE){
    eW <- eigen(W, symmetric=symmetric)
    d <- eW$values
    if (min(d) <0  && abs(min(d))>sqrt(.Machine$double.eps))
      stop("'W' is not positive definite")
    else d[d<=0]<- ifelse(inverse, Inf, 0)
    A <- diag(d^ifelse(inverse, -0.5, 0.5)) %*% t(eW$vector)
    A # t(A)%*%A = W^{-1}
  }
  
  variance.decomposition.old <- function(y, covar, G) {
    EE <- diag(nrow(G))
    if (is.null(covar)) {
      vc <- estVC(y=pheno, v=list(AA=G, DD=NULL, HH=NULL, AD=NULL, MH=NULL, EE=EE))
    } else {
      vc <- estVC(y=pheno, covar, v=list(AA=G, DD=NULL, HH=NULL, AD=NULL, MH=NULL, EE=EE))
    }  
    V <- vc$par["AA"]*G + vc$par["EE"]*EE
    A <- W.inv(V)
    return(list(A=A, V=V, var.g = vc$par["AA"], var.e = vc$par["EE"]))
  }
  
  variance.decomposition <- function(y, covar, G) {
    fit <- regress(y~.,~G, data=as.data.frame(cbind(y=y, covar))) 
    V <- fit$sigma["G"] * G + fit$sigma["In"] * diag(nrow(G))
    A <- W.inv(V)
    return(list(A=A, V=V, var.g = fit$sigma["G"], var.e = fit$sigma["In"]))
  }
  
  nonNA <- !is.na(pheno)
  
  procedure <- match.arg(procedure)
  
  # linear model
  if (procedure == "LM") G <- A <- NULL
  
  # linear mixed model
  if (procedure == "LMM") {
    if (is.null(G)) {
      warning("Genetic similarity matrix G should be specified, use 'gensim.matrix' function.")
      G <- gensim.matrix(geno, ...)
    }  
    A <- variance.decomposition(pheno[nonNA], covar[nonNA,], G[nonNA,nonNA])$A
  }  
  
  # linear mixed model, leave one out
  if (procedure == "LMM-L1O") {
    if (is.null(geno$chromosomes) | is.null(geno$markers$chr) | is.null(geno$markers$pos))
      stop("For LMM-L1O procedure, markers MUST be mapped to chromosomes.")
    if (is.null(G)) {
      warning("Genetic similarity matrix G should be specified, use 'gensim.matrix' function.")
      G <- list()
      for (c in geno$chromosomes$chr)
        G[[c]] <- gensim.matrix(geno, skip.chr=c, ...)
    }
    
    A <- list()
    for (c in geno$chromosomes$chr) {
      A[[c]] <- variance.decomposition(pheno[nonNA], covar[nonNA,], G[[c]][nonNA,nonNA])$A
    }
  } 
  
  nind <- sum(nonNA) # number of individuals
  output <- data.frame(chr=geno$markers$chr, pos=geno$markers$pos, lod=rep(0, nrow(geno$markers)), row.names=geno$markers$marker)
  
  for (c in geno$chromosomes$chr) {
    if (is.list(A)) A.chr <- A[[c]] else A.chr <- A
    
    # if A.chr is given, transform variables
    if (!is.null(A.chr)) ynew <- A.chr %*% pheno[nonNA] else ynew <- pheno[nonNA]
    if (!is.null(A.chr)) covarnew <- A.chr %*% cbind(rep(1,nind),covar[nonNA,]) else covarnew <- cbind(rep(1,nind),covar[nonNA,])
    
    # null model
    rss0 <- sum(lsfit(y=ynew,x=covarnew,intercept=FALSE)$residuals^2)
    
    for (i in which(geno$markers$chr == c)){
      if (!is.null(A.chr)) xnew <- A.chr %*% cbind(geno$probs[nonNA,,i],covar[nonNA,]) else xnew <- cbind(geno$probs[nonNA,,i], covar[nonNA,])
      rss1 <- sum(lsfit(y=ynew,x=xnew,intercept=FALSE)$residuals^2)
      output$lod[i] <- nind/2 * (log10(rss0) - log10(rss1))
    } 
  }
  
  class(output) <- c("scanone", "data.frame")
  return(output) 
  
}