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

  \item{package}{/ method to do variant component
  calculations}
}
\value{
4-elements list
}
\description{
Decompose variance into genetic part and a random noise
}
\details{
Currently, \code{regress} function is much faster than
\code{QTLRel::estVC} but it prints
}
\examples{
cross <- sim.cross.geno(250, nmar=10)
cross$pheno <- sim.cross.pheno(0.5, cross)
heritability(cross)
}
\keyword{manip}

