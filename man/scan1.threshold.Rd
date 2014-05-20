\name{scan1.threshold}
\alias{scan1.threshold}
\title{Significance Threshold of Genome Scan}
\usage{
scan1.threshold(geno, pheno, pheno.cols = 1, covar = NULL,
  procedure = c("LM", "LMM", "LMM-L1O"), G, n.perm = 100, alpha = 0.05,
  subjects = seq(geno$subjects), markers = seq(NROW(geno$markers)), ...)
}
\arguments{
  \item{geno}{genotype probabilities
  ("\code{genotype.probs}" or "\code{cross}" object)}

  \item{pheno}{data frame with phenotypes}

  \item{pheno.cols}{selection of phenotype's column(s)}

  \item{covar}{(additive) covariates}

  \item{procedure}{procedure to do the inference, see
  Details}

  \item{G}{genetic similarity matrix or a list of genetic
  similarity matrices or \code{NULL}}

  \item{subjects}{subseting of subjects}

  \item{markers}{subseting of markers}

  \item{n.perm}{number of permutations}

  \item{alpha}{level of significance}

  \item{...}{parameters passed to \code{genrel.matrix}}
}
\value{
numeric
}
\description{
Significance Threshold of Genome Scan
}
\details{
Currently, three procedures are implemented: linear model
(LM), linear mixed model (LMM) and linear mixed model -
leave one chromosome out (LMM-L1O).
}
\examples{
data(fake.f2, package="qtl")
fake.f2 <- calc.genoprob(fake.f2)

# warning, n.perm=10 is too low for practical purposes (but fast)
scan1.threshold(fake.f2, procedure="LM", n.perm=10)
scan1.threshold(fake.f2, procedure="LMM", n.perm=10)
scan1.threshold(fake.f2, procedure="LMM-L1O", n.perm=10)
}
\keyword{manip}

