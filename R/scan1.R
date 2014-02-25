#' Genome Scan With a Single QTL Model
#' 
#' @param cross "\code{cross}" object
#' @param pheno.cols selection of phenotype's column(s)
#' @param procedure procedure to do the inference
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

scan1 <- function(cross, pheno.cols=1, procedure="scanone", ...) {
  
  cross <- calc.genoprob(cross)
  output <- qtl::scanone(cross, method="hk", pheno.col=pheno.cols)
  
  if (procedure == "scanone") {
    return(output)
  }
  
  G <- genrel.matrix(cross, ...)   
  
  if (procedure == "scanOne" || procedure == "regress-qtl") {
    for (p in 1:length(pheno.cols)) {
      Y <- cross$pheno[,pheno.cols[p]]
      EE <- diag(nind(cross))
      vc <- estVC(y=Y, v=list(AA=G, DD=NULL, HH=NULL, AD=NULL, MH=NULL, EE=EE))
      if (procedure == "scanOne") {
        output[,p+2] <- scanOne(y=Y, gdat=extract.geno(cross), vc=vc)$p / (2*log(10))
      } else {
        # variance matrix
        V <- vc$par[2]*G + vc$par[3]*EE
        
        # calculate A = V^-(1/2)
        ei <- eigen(solve(V))
        A <- ei$vectors %*% diag(sqrt(ei$values)) %*% t(ei$vectors)
        
        # Y is tranformed to remove a correlation structure
        Y.transformed <- A %*% (Y - vc$par[1])
        cross$pheno[,pheno.cols[p]] <- Y.transformed
        output[,p+2] <- scanone(cross, method="hk", pheno.col = pheno.cols[p])[,3]
      }  
    }
    return(output)
  }
  
  
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
}