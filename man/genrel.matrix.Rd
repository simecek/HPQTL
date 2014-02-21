\name{genrel.matrix}
\alias{genrel.matrix}
\title{Genetic Relationship Matrix}
\usage{
genrel.matrix(geno, method = "additive", expected.means = NULL)
}
\arguments{
  \item{geno}{\code{cross} object or a genotype matrix
  (\code{nind} rows, \code{totalmar} columns)}

  \item{method}{how to calculate similarity matrix
  ('additive', 'additive-diag1', 'kinship')}

  \item{expected.means}{a vector of expected marker means,
  if \code{NULL} means are estimated from data}
}
\value{
square matrix (\code{nind} rows, \code{nind} columns)
}
\description{
Genetic Relationship Matrix
}
\examples{
cross <- sim.cross.geno(250, nmar=10)
geno <- extract.geno(cross)
genrel.matrix(geno)
genrel.matrix(geno, expected.means = c(rep(2,190), rep(1.5,10)))
}
\keyword{manip}

