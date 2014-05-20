\name{heritability}
\alias{heritability}
\title{Estimate Heritability}
\usage{
heritability(geno, pheno, pheno.cols = 1, covar = NULL, se = FALSE,
  G = NULL, ...)
}
\arguments{
  \item{geno}{genotype probabilities
  ("\code{genotype.probs}" or "\code{cross}" object)}

  \item{pheno}{data frame with phenotypes}

  \item{pheno.cols}{selection of phenotype's column(s)}

  \item{covar}{(additive) covariates}

  \item{G}{genetic similarity matrix or a list of genetic
  similarity matrices or \code{NULL}}

  \item{se}{if se=TRUE standard error is estimates}

  \item{...}{parameters passed to \code{gensim.matrix}}
}
\value{
numeric
}
\description{
Estimate Heritability
}
\keyword{manip}

