\name{extract.geno}
\alias{extract.geno}
\title{Extract Genotype Probabilities}
\usage{
extract.geno(cross)
}
\arguments{
  \item{cross}{"\code{cross}" object}
}
\value{
\code{genotype.probs} is a list of three components
\itemize{ \item \code{probs} \item \code{calls} \item
\code{markers} }
}
\description{
Extract Genotype Probabilities
}
\details{
Extract genotype prob. array from "\code{qtl::cross}"
object. Stored as a list
}
\examples{
cross <- sim.cross.geno(250, nmar=10)
extract.geno(cross)
}
\keyword{manip}

