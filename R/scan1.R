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

scan1 <- function(pheno, geno, covar=NULL, procedure=c("lm","lmmi","lmme","user-defined"), G=NULL, A=NULL, ...) {

  W.inv<- function(W, symmetric=TRUE,inverse=TRUE){
    eW <- eigen(W, symmetric=symmetric)
    d <- eW$values
    if (min(d) <0  && abs(min(d))>sqrt(.Machine$double.eps))
      stop("'W' is not positive definite")
    else d[d<=0]<- ifelse(inverse, Inf, 0)
    A <- diag(d^ifelse(inverse, -0.5, 0.5)) %*% t(eW$vector)
    A # t(A)%*%A = W^{-1}
  }
  
  procedure <- match.arg(procedure)
  
  if (procedure == "lm")   G <- NULL
  if (procedure == "lmmi") G <- genrel.matrix(geno, ...) 
  if (procedure == "lmme") {
    G <- list()
    for (i in 1:length(geno$probs)) {
      G[[i]] <- genrel.matrix(geno, skip.chr=i, ...)
    }
  } 
  
  if (is.null(A)) {
    if (is.null(G)) A <- NULL
    if (is.matrix(G)) {
      EE <- diag(nrow(G))
      if (is.null(covar)) {
        vc <- estVC(y=pheno[,1], v=list(AA=G, DD=NULL, HH=NULL, AD=NULL, MH=NULL, EE=EE))
      } else {
        vc <- estVC(y=pheno[,1], covar, v=list(AA=G, DD=NULL, HH=NULL, AD=NULL, MH=NULL, EE=EE))
      }  
      V <- vc$par["AA"]*G + vc$par["EE"]*EE
      A <- W.inv(V)
    }
    if (is.list(G)) {
      A <- list()
      for (i in 1:length(geno$probs)) {
        EE <- diag(nrow(G[[i]]))
        if (is.null(covar)) {
          vc <- estVC(y=pheno[,1], v=list(AA=G[[i]], DD=NULL, HH=NULL, AD=NULL, MH=NULL, EE=EE))
        } else {
          vc <- estVC(y=pheno[,1], covar, v=list(AA=G[[i]], DD=NULL, HH=NULL, AD=NULL, MH=NULL, EE=EE))
        }  
        V <- vc$par["AA"]*G[[i]] + vc$par["EE"]*EE
        A[[i]] <- W.inv(V)
      }  
    }
  }  

  nind <- nrow(geno$probs[[1]]) # number of individuals
  output <- NULL
  
  for (i in 1:length(geno$probs)) {
    if (is.list(A)) A.chr <- A[[i]] else A.chr <- A
    
    if (!is.null(A.chr)) ynew <- A.chr %*% pheno else ynew <- pheno
    if (!is.null(A.chr)) covarnew <- A.chr %*% cbind(rep(1,nind),covar) else covarnew <- cbind(rep(1,nind),covar)
    rss0 <- sum(lsfit(y=ynew,x=covarnew,intercept=FALSE)$residuals^2)
    
    rss1 <- rep(0.0, nrow(geno$marker[[i]]))
    for (j in 1:nrow(geno$marker[[i]])){
      if (!is.null(A.chr)) xnew <- A.chr %*% cbind(geno$probs[[i]][,,j],covar) else xnew <- cbind(geno$probs[[i]][,,j],covar)
      rss1[j] <- sum(lsfit(y=ynew,x=xnew,intercept=FALSE)$residuals^2)
    }
    
    output.chr <- geno$markers[[i]][,-1]
    rownames(output.chr) <- geno$markers[[i]][,1]
    output.chr$lod <- nind/2 * (log10(rss0) - log10(rss1))
    output <- rbind(output, output.chr)
    class(output) <- c("scanone", "data.frame")
  }
  return(output)
  
  if (procedure == "scanOne-per-chr") {
    #extract genotype
    geno <- extract.geno(cross)
    
    # get a vector of markers' chromosomes
    nmar <- nmar(cross)
    chrs <- names(nmar)
    marchrs <- c()
    for (i in 1:length(nmar))
    marchrs <- c(marchrs, rep(names(nmar)[i], nmar[i]))
    
    for (p in 1:length(pheno.cols)) {
      Y <- cross$pheno[,pheno.cols[p]]
      EE <- diag(nind(cross))
           
      for (c in chrs) {
        G.minus.c <- genrel.matrix(geno[,marchrs!=c])
        vc.c <- estVC(y=Y, v=list(AA=G.minus.c, DD=NULL, HH=NULL, AD=NULL, MH=NULL, EE=EE))
        output[marchrs==c,p+2] <- scanOne(y=Y, gdat=geno[,marchrs==c], vc=vc.c)$p / (2*log(10))
      }
    }
    
    return(output)
  }

  if (procedure == "blup") {
    Z <- extract.geno(cross)
    W <- Z %*% t(Z)
    E <- diag(nrow(W))
    
    for (p in 1:length(pheno.cols)) {
      Y <- cross$pheno[,pheno.cols[p]]
      fit <- regress(Y~1, ~W)
      Sigma <- fit$sigma["W"] * W + fit$sigma["In"] * E
      u.est <- sqrt(fit$sigma["W"]) * t(Z) %*% (solve(Sigma, Y - fit$fitted)) 
      output[,p+2] <- abs(u.est - median(u.est))
    }
    return(output)
  }
  
  # not working
  if (procedure == "vc-qtl-per-chr") {
    
    #extract genotype
    geno <- extract.geno(cross)
    
    # get a vector of markers' chromosomes
    nmar <- nmar(cross)
    chrs <- names(nmar)
    marchrs <- c()
    for (i in 1:length(nmar))
      marchrs <- c(marchrs, rep(names(nmar)[i], nmar[i]))
    
    for (p in 1:length(pheno.cols)) {
      Y <- cross$pheno[,pheno.cols[p]]
      EE <- diag(nind(cross))
      
      for (c in chrs) {
        G.minus.c <- genrel.matrix(geno[,marchrs!=c])
        vc <- estVC(y=Y, v=list(AA=G.minus.c, DD=NULL, HH=NULL, AD=NULL, MH=NULL, EE=EE))
        V.c <- vc$par[2]*G.minus.c + vc$par[3]*EE
        
        # calculate A = V^-(1/2)
        ei <- eigen(solve(V.c))
        A <- ei$vectors %*% diag(sqrt(ei$values)) %*% t(ei$vectors)
        
        # use scanone to transformed Y
        Y.transformed <- A %*% (Y - vc$par[1])
        cross$pheno[,p] <- Y.transformed
        output[marchrs==c,p+2] <- scanone(cross, method="hk", pheno.col = pheno.cols[p])[marchrs==c,3]
      }
    }  
    return(output)  
  }
  
  
}