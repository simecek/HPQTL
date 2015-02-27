#' Import qtl / DOQTL cross
#'
#' @details Extract genotype prob. array from "\code{qtl::cross}" object. 
#'   
#' @param cross "\code{cross}" object or DO qtl 
#' @param check.output check consistency of results
#'
#' @return \code{genotype.probs} object is a list of three components 
#' \itemize{ 
#'   \item \code{probs} 3-dimensional array of genotype probabilities (subjects x calls x markers) or a list of such arrays
#'   \item \code{calls} possible genotype calls 
#'   \item \code{markers}
#' } 
#' 
#' @keywords manip
#'
#' @export 
#'
#' @return \code{geno} object
#'   
#' @examples
#' data(fake.f2, package="qtl")
#' fake.f2 <- calc.genoprob(fake.f2)
#' geno <- extract.geno(fake.f2)

extract.geno <- function(cross, check.output = TRUE) {

  # currently only cross (qtl package) or DO genotype probs implemented
  stopifnot("cross" %in% class(cross) | "array" %in% class(cross))
  
  if ("cross" %in% class(cross)) {
  
    # check that cross is either bc or f2
    if (class(cross)[1] != "bc" & class(cross)[1] != "f2")
      warning("The class of cross should be either 'bc' or 'f2'.")
      
    # check if calc.genoprob has been called
    if (!("prob" %in% names(cross$geno[[1]]))) {
      warning("First running calc.genoprob.")
      cross <- calc.genoprob(cross)
    }

    # chromosome types ("A" or "X")
    chrtype <- sapply(cross$geno, class)
    
    # for chrX data, revise probs
    if(any(chrtype =="X")) {
      for(i in which(chrtype=="X"))
        cross$geno[[i]]$prob <- qtl:::reviseXdata(class(cross)[1], "simple", cross$pheno, prob=cross$geno[[i]]$prob)
    }
    
    # extract genotype probabilities 
    list.of.probs <- lapply(cross$geno, function(x) aperm(x$prob, c(1,3,2))  )
    
    # guess if first column in phenotype table could be used as subject id
    if (class(cross$pheno[,1]) == "character" & length(unique(cross$pheno[,1])) == nrow(cross$pheno)) {
      subjects = cross$pheno[,1]
    } else {
      subjects = as.character(1:nrow(cross$pheno))
    }
    
    # table of markers and their positions
    markers <- data.frame(marker = do.call("c", lapply(cross$geno, function(x) names(x$map))),
                          chr = do.call("c", mapply(rep, x = names(nmar(cross)), each = nmar(cross), SIMPLIFY = FALSE)),
                          pos = do.call("c", lapply(cross$geno, function(x) as.vector(x$map))), 
                          stringsAsFactors = FALSE)
    
    # if number of call differes
    if (length(unique(sapply(list.of.probs, ncol)))!=1) {
      warning("Number of calls differ across chromosomes.")
      max.calls <- max(sapply(list.of.probs, ncol))
      for (i in which(sapply(list.of.probs, ncol) < max.calls)) {
        miss.calls <- max.calls - dim(list.of.probs[[i]])[2] 
        zeros <- array(rep(0,dim(list.of.probs[[i]])[1] * miss.calls * dim(list.of.probs[[i]])[3]))
        dim(zeros) <- c(dim(list.of.probs[[i]])[1], miss.calls,  dim(list.of.probs[[i]])[3])
        list.of.probs[[i]] <- abind(list.of.probs[[i]], zeros, along = 2)
      }
    }
    
    # see help for 'geno' class
    geno <- list(probs = abind(list.of.probs, along=3), 
                 subjects = subjects,    
                 calls = dimnames(cross$geno[[1]]$prob)[[3]],
                 markers = markers,
                 chromosomes = data.frame(chr = names(cross$geno), type = as.vector(chrtype), stringsAsFactors = FALSE))
  }
  
  # 3-dim DO array with attribute markers
  if ("array" %in% class(cross)) {
    
    subjects <- dimnames(cross)[[1]]
    calls <- dimnames(cross)[[2]]
    markers <- attr(cross, "markers")
    chromosomes <- attr(cross, "chromosomes")
    
    # if markers or chromosomes are missing, try to get provisional ids
    if (is.null(markers)) {
      warning("Attribute 'markers' is missing")
      markers <- data.frame(marker = dimnames(cross)[[3]], stringsAsFactors = FALSE)
      
    }
    if (is.null(chromosomes) & !is.null(markers$chr)) {
      chrs <- unique(markers$chr)
      chromosomes <- data.frame(chr=chrs, type=ifelse(chrs=="X","X","A"), stringsAsFactors = FALSE)
    }

    # remove attributes
    #attributes(cross) <- NULL
    
    geno <- list(probs = cross, 
                    subjects = subjects,    
                    calls = calls,
                    markers = markers,
                    chromosomes = chromosomes)
  }
    
  class(geno) <- c("genotype.probs", "list")
  if (check.output) check.genotype.probs(geno)
  geno
}