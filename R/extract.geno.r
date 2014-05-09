#' Import qtl/DOQTL cross
#'
#' @details Extract genotype prob. array from "\code{qtl::cross}" object. 
#'   
#' @param cross "\code{cross}" object
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
#' cross <- sim.cross.geno(250, nmar=10)
#' extract.geno(cross)

extract.geno <- function(cross) {

  # currently only cross (qtl package) or DO genotype probs implemented
  stopifnot("cross" %in% class(cross) | "array" %in% class(cross))
  
  if ("cross" %in% class(cross)) {
  
    # check that calc.genoprob has been called
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
                                  chr = do.call("c", mapply(rep, x = names(nmar(cross)), each = nmar(cross))),
                                  pos = do.call("c", lapply(cross$geno, function(x) as.vector(x$map))))
    
    # see help for 'geno' class
    geno <- list(probs = abind(list.of.probs, along=3), 
                 subjects = subjects,    
                 calls = dimnames(cross$geno[[1]]$prob)[[3]],
                 markers = markers,
                 chromosomes = data.frame(chr = names(cross$geno), type = as.vector(chrtype)))
  }
  
  # DO array with 
  if ("array" %in% class(cross)) {
    
    subjects <- dimnames(cross)[[1]]
    calls <- dimnames(cross)[[2]]
    markers <- attr(cross, "markers")
    chromosomes <- attr(cross, "chromosomes")
    
    # if markers or chromosomes are missing, try to get provisional
    if (is.null(markers)) {
      warning("Attribute 'markers' is missing")
      markers <- data.frame(marker = dimnames(cross)[[3]])
      
    }
    if (is.null(chromosomes) & !is.null(markers$chr)) {
      chrs <- unique(markers$chr)
      chromosomes <- data.frame(chr=chrs, type=ifelse(chrs=="X","X","A"))
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
  #check.genotype.probs(geno)
  geno
}