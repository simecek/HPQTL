#' Decompose variance into genetic part and a random noise
#' 
#' @param y a vector with phenotype
#' @param covar a matrix with covariates
#' @param G genetic similarity matrix
#' @param package / method to do variant component calculations
#' 
#' @details Currently, \code{regress} function is much faster than \code{QTLRel::estVC}
#' but it prints 
#' 
#' 
#' @return 4-elements list
#' 
#' @keywords manip
#'
#' @examples
#' cross <- sim.cross.geno(250, nmar=10)
#' cross$pheno <- sim.cross.pheno(0.5, cross)
#' heritability(cross)

variance.decomposition <- function(y, covar, G, package = c("regress", "QTLRel"), ...) {
  package = match.arg(package)
  
  if (package=="regress") {
    fit <- regress(y~.,~G, data=as.data.frame(cbind(y=y, covar))) 
    V <- fit$sigma["G"] * G + fit$sigma["In"] * diag(nrow(G))
    A <- half.inv(V)
    return(list(V=V, A=A, var.g = fit$sigma["G"], var.e = fit$sigma["In"]))
  } else {
    EE <- diag(nrow(G))
    if (is.null(covar)) {
      vc <- estVC(y=y, v=list(AA=G, DD=NULL, HH=NULL, AD=NULL, MH=NULL, EE=EE))
    } else {
      vc <- estVC(y=y, covar, v=list(AA=G, DD=NULL, HH=NULL, AD=NULL, MH=NULL, EE=EE))
    }  
    V <- vc$par["AA"]*G + vc$par["EE"]*EE
    A <- half.inv(V)
    return(list(V=V, A=A, var.g = vc$par["AA"], var.e = vc$par["EE"]))
  }  
}