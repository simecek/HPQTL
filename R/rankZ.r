#' RankZ Normalization transform
#'
#' Return quantiles normalized (rank-Z) transformed vector.  
#' 
#' @param x numeric vector
#'
#' @keywords manip
#' 
#' @return Normalized vector. 
#'
#' @examples
#'
#' x <- 2 ^ rnorm(1000) 
#' hist(x)
#' 
#' x_rankz <- rankZ(x)
#' hist(x_rankz)

rankZ <- function(x) {
  y <- rank(x)
  y[is.na(x)] <- NA
  qnorm(y / (sum(!is.na(x))+1))
}