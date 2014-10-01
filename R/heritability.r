#' Estimate Heritability
#' 
#' @param y a vector with phenotype
#' @param covar a matrix with covariates, Intercept should be included
#' @param G genetic similarity matrix
#' @param se if TRUE standard error is estimated
#' @param var.decompose if TRUE total variability is decomposed as into three parts: fixed effects, animal effect and noise
#' @param ... parameters passed to \code{gensim.matrix}
#' @return numeric
#' 
#' @keywords manip
#'
#' @export

heritability <- function(y, covar, G, se=FALSE, var.decompose = FALSE, ...) {
  
  # if 'covar' is missing, make it intercept
  if (missing(covar))
    covar <- rep(1, nrow(G))
  
  # if 'G' is missing then it should be estimated from genotype  
  if (missing(G)) {
    G <- gensim.matrix(geno, ...)
  }
    
  # fit the mixed model
  rg.fit <- regress(y~covar, ~G, pos=TRUE)
  
  # estimate heritability
  h2 <- as.numeric(rg.fit$sigma[1] / sum(rg.fit$sigma))

  # calculate heritability se
  if (se) {
    v <- c(rg.fit$sigma[2], -rg.fit$sigma[1]) / sum(rg.fit$sigma^2)
    h2.se <- rbind(v) %*% rg.fit$sigma.cov %*% cbind(v)
    attr(h2, "se") <- as.numeric(sqrt(h2.se))
  }  
  
  # if total decomposition is asked, estimate it
  if (var.decompose) {
    X <- model.matrix(~covar)
    # drop Intercept if not estimated
    if (ncol(X) == length(rg.fit$beta) + 1) {
      X <- X[,-1]
    }
    var.fixed <- var(X %*% cbind(rg.fit$beta))
    var.total <- var.fixed + sum(rg.fit$sigma)
    h2.var.decompose <- c(var.fixed, rg.fit$sigma) / var.total
    names(h2.var.decompose) <- c("fixed","background","noise")
    attr(h2, "var.decompose") <- h2.var.decompose
  }
  
  return(h2)
}