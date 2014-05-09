\name{import.cross}
\alias{import.cross}
\title{Import qtl/DOQTL cross}
\usage{
import.cross(cross)
}
\arguments{
  \item{cross}{"\code{cross}" object}
}
\value{
\code{genotype.probs} object is a list of three components
\itemize{ \item \code{probs} 3-dimensional array of
genotype probabilities (subjects x calls x markers) or a
list of such arrays \item \code{calls} possible genotype
calls \item \code{markers} }

\code{hpcross} object
}
\description{
Import qtl/DOQTL cross
}
\details{
Extract genotype prob. array from "\code{qtl::cross}"
object.
}
\examples{
cross <- sim.cross.geno(250, nmar=10)
extract.geno(cross)
}
\keyword{manip}

