\name{heritability}
\alias{heritability}
\title{Estimate Heritability}
\usage{
heritability(cross, pheno.cols = 1, se = FALSE, ...)
}
\arguments{
  \item{cross}{"\code{cross}" object}

  \item{pheno.cols}{selection of phenotype's column(s)}

  \item{se}{if se=TRUE standard error is estimates}

  \item{...}{parameters passed to \code{genrel.matrix}}
}
\value{
numeric
}
\description{
Estimate Heritability
}
\examples{
cross <- sim.cross.geno(250, nmar=10)
cross$pheno <- sim.cross.pheno(0.5, cross)
heritability(cross, se=TRUE)
}
\keyword{manip}

