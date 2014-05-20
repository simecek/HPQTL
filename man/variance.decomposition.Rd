\name{variance.decomposition}
\alias{variance.decomposition}
\title{Decompose variance into genetic part and a random noise}
\usage{
variance.decomposition(y, covar, G, package = c("regress", "QTLRel"), ...)
}
\arguments{
  \item{y}{a vector with phenotype}

  \item{covar}{a matrix with covariates}

  \item{G}{genetic similarity matrix}

  \item{package}{package to do variant component
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
\code{QTLRel::estVC}.
}
\examples{
data(fake.f2, package="qtl")
fake.f2 <- calc.genoprob(fake.f2)

geno <- extract.geno(fake.f2)
G <- gensim.matrix(geno)
HPQTL:::variance.decomposition(fake.f2$pheno[,1], NULL, G)
}
\keyword{manip}

