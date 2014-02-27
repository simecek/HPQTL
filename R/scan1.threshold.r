#' Significance Threshold of Genome Scan
#'  
#' @param cross "\code{cross}" object
#' @param pheno.cols selection of phenotype's column(s)
#' @param procedure procedure to do the inference
#' @param n.perm number of permutations
#' @param alpha level of significance 
#' @param ... parameters passed to \code{genrel.matrix}
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
#' heritability(cross)

scan1.threshold <- function(cross, pheno.cols=1, procedure="scanone", n.perm=100, alpha=0.05, ...) {
  
  cross <- calc.genoprob(cross)
  
  if (procedure == "scanone") {
    trhold <- summary(scanone(cross, n.perm=n.perm, pheno.col=pheno.cols, method="hk", verbose=FALSE), alpha=alpha)
    output <- as.numeric(trhold)
  }
    
  if (procedure == "scanOne") {
    
    geno <- extract.geno(cross)
    G <- genrel.matrix(geno, ...)
    EE <- diag(nind(cross))
    vc <- list()
    for (p in 1:length(pheno.cols))
      Y <- cross$pheno[,pheno.cols[p]]
      vc[[p]] <- estVC(y=Y, v=list(AA=G, DD=NULL, HH=NULL, AD=NULL, MH=NULL, EE=EE))
    maxs <- matrix(0, n.perm, length(pheno.cols))
    
    for (i in 1:n.perm) {
      pheno.perm <- sample(nind(cross))
      for (p in 1:length(pheno.cols)) {
        Y <- cross$pheno[pheno.perm,pheno.cols[p]]
        maxs[i,p] <- max(scanOne(y=Y, gdat=geno, vc=vc[[p]])$p) / (2*log(10))
      }  
    }
    
    output <- apply(maxs, 2, quantile, probs = 1 - alpha)
  }
  
  if (procedure == "scanOne-per-chr") {
    #extract genotype
    geno <- extract.geno(cross)
    
    # get a vector of markers' chromosomes
    nmar <- nmar(cross)
    chrs <- names(nmar)
    marchrs <- c() # chromosome of markers 
    for (i in 1:length(nmar))
      marchrs <- c(marchrs, rep(names(nmar)[i], nmar[i]))
    output <- c()
    
    for (p in 1:length(pheno.cols)) {
      Y <- cross$pheno[,pheno.cols[p]]
      EE <- diag(nind(cross))
      Intercept <- cbind(rep(1,nind(cross)))
      perms <- matrix(nrow = length(chrs),ncol = n.perm)
      
      for (c in chrs) {
        G.minus.c <- genrel.matrix(geno[,marchrs!=c])
        vc.c <- estVC(y=Y, v=list(AA=G.minus.c, DD=NULL, HH=NULL, AD=NULL, MH=NULL, EE=EE))
        V <- vc.c$par[2]*G.minus.c + vc.c$par[3]*EE
        Y.perm <- mvnpermute(Y, Intercept, V, n.perm)
        for (i in 1:n.perm) {
          lod <- scanOne(y=Y.perm[,i], gdat=geno[,marchrs==c], vc=V)$p / (2*log(10))
          perms[which(chrs==c),i] <- max(lod)
        }
      }
      
      trhold <- quantile(apply(perms, 2, max), 1-alpha)
      output <- c(output, trhold)
    }

  }  
  
  if (procedure == "vc-qtl") {
    geno <- extract.geno(cross)
    G <- genrel.matrix(geno, ...)
    EE <- diag(nind(cross))
    
    for (p in 1:length(pheno.cols)) {
      Y <- cross$pheno[,pheno.cols[p]]
      vc <- estVC(y=Y, v=list(AA=G, DD=NULL, HH=NULL, AD=NULL, MH=NULL, EE=EE))
      V <- vc$par[2]*G + vc$par[3]*EE
      ei <- eigen(solve(V))
      A <- ei$vectors %*% diag(sqrt(ei$values)) %*% t(ei$vectors)
      Y.transformed <- A %*% (Y - vc$par[1])
      cross$pheno[,pheno.cols[p]] <- Y.transformed
    }
    
    trhold <- summary(scanone(cross.observed, n.perm=n.perm, pheno.col=pheno.cols, method="hk", verbose=FALSE), alpha=alpha)
    output <- as.numeric(trhold)
  }
  
  return(as.numeric(output))
}