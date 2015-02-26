#' Genetic Similarity Matrix
#' 
#'   
#' @param geno genotype probabilities ("\code{genotype.probs}" or "\code{cross}" object)
#' @param procedure procedure of scan1 the G is calculated for
#' @param subjects subseting of subjects
#' @param markers subseting of markers
#' @param gensim.normalization method to be used for genetic similarity matrix normalization   
#'
#' @return square matrix (\code{nind} rows, \code{nind} columns)
#' 
#' @keywords manip
#'
#' @export
#' 
#' @examples
#' data(fake.f2, package="qtl")
#' fake.f2 <- calc.genoprob(fake.f2)
#' 
#' geno <- extract.geno(fake.f2)
#' G <- gensim.matrix(geno)
#' Glist <- gensim.matrix(geno, procedure="LOCO")

gensim.matrix <- function(geno, procedure = c("LMM", "LOCO", "POMOLOCO", "LM"),  
                          gensim.normalization="sample-variance", verbose = FALSE) {
  
  # match arguments
  procedure = match.arg(procedure)
  
  # if not genotype.probs then extract genotype probs
  if (!("genotype.probs" %in% class(geno))) geno <- extract.geno(geno)
  check.genotype.probs(geno)
  
  # for linear model no genetic similarity matrix is needed
  if (procedure == "LM") return(NULL)

  # calculate GSM per chromosome
  GSM.per.chr <- list()
  for (c in geno$chromosomes$chr) {
    if (verbose) message(paste("Processing GSM of "))
    tmp <- geno$probs[,,geno$markers$chr==c]
    dim(tmp) <- c(dim(tmp)[1], dim(tmp)[2] * dim(tmp)[3])
    GSM.per.chr[[c]] <- tcrossprod(tmp)
  }
  
  if (procedure == "LMM") {
    G = Reduce("+", GSM.per.chr)
    return(normalize.matrix(G, method = gensim.normalization))
  }
  
  if (procedure == "LOCO") {
    Glist <- list()
    for (c in geno$chromosomes$chr)
      Glist[[c]] <- Reduce("+", GSM.per.chr[geno$chromosomes$chr!=c])
    return(normalize.matrix(Glist, method = gensim.normalization))
  }
  
  if (procedure == "POMOLOCO") {
    return(GSM.per.chr)
  }  
    
  stop("Something wrong happened.")
}