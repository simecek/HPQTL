#' Decompose variance into genetic part and a random noise
#' 
#' @param y a vector with phenotype
#' @param covar a matrix with covariates
#' @param G genetic similarity matrix
#' @param package package to do variant component calculations
#' 
#' @details Currently, \code{regress} function is much faster than \code{QTLRel::estVC}.
#'  
#' @return 4-elements list
#' 
#' @keywords manip
#'
#' @examples
#' data(fake.f2, package="qtl")
#' fake.f2 <- calc.genoprob(fake.f2)
#' 
#' geno <- extract.geno(fake.f2)
#' G <- gensim.matrix(geno)
#' variance.decomposition(fake.f2$pheno[,1], NULL, G)
#' 
#' @export


variance.decomposition <- function(y, covar, G, package = c("regress", "QTLRel"), ...) {
  package = match.arg(package)
  
  if (missing(covar)) covar <- cbind(rep(1, nrow(G)))
  
  if (package=="regress") {
    fit <- regress(y~covar, ~G, pos=c(TRUE, TRUE)) 
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
