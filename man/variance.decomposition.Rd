\name{variance.decomposition}
\alias{variance.decomposition}
\title{Decompose variance into genetic part and a random noise}
\usage{
variance.decomposition(y, covar, G, package = c("regress", "QTLRel"))
}
\arguments{
  \item{y}{a vector with phenotype}

  \item{covar}{a matrix with covariates}

  \item{G}{genetic similarity matrix}
}
\value{
4-elements list
}
\description{
Decompose variance into genetic part and a random noise
}
\examples{
cross <- sim.cross.geno(250, nmar=10)
cross$pheno <- sim.cross.pheno(0.5, cross)
heritability(cross)
}
\keyword{manip}

