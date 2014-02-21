\name{sim.cross.pheno}
\alias{sim.cross.pheno}
\title{Generate Random Phenotype with Given Heritability}
\usage{
sim.cross.pheno(h2, cross, background = "GRM", ...)
}
\arguments{
  \item{h2}{heritability}

  \item{cross}{\code{cross} object}

  \item{background}{method to generate genetic background}

  \item{...}{parameters passed to \code{genrel.matrix}}
}
\value{
"\code{cross}" object
}
\description{
Generate Random Phenotype with Given Heritability
}
\examples{
cross <- sim.cross.geno(250, nmar=10)
cross$pheno <- sim.cross.pheno(0.5, cross)
}
\keyword{manip}

