\name{extract.geno}
\alias{extract.geno}
\title{Extract Genotype Matrix}
\usage{
extract.geno(cross)
}
\arguments{
  \item{cross}{"\code{cross}" object}
}
\value{
matrix (\code{nind} rows, \code{totalmar} columns)
}
\description{
Extract Genotype Matrix
}
\examples{
cross <- sim.cross.geno(250, nmar=10)
extract.geno(cross)
}
\keyword{manip}

