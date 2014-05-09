\name{check.genotype.probs}
\alias{check.genotype.probs}
\title{Check \code{genotype.probs} object}
\usage{
check.genotype.probs(geno)
}
\arguments{
  \item{geno}{a \code{genotype.probs} object}
}
\value{
\code{NULL}
}
\description{
Check \code{genotype.probs} object
}
\examples{
cross <- sim.cross.geno(250, nmar=10)
geno <- extract.geno(cross)
check.genotype.probs(geno)
}
\keyword{manip}

