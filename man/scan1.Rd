\name{scan1}
\alias{scan1}
\title{Genome Scan With a Single QTL Model}
\usage{
scan1(cross, pheno.cols = 1, procedure = "scanone", G = NULL, ...)
}
\arguments{
  \item{cross}{"\code{cross}" object}

  \item{pheno.cols}{selection of phenotype's column(s)}

  \item{procedure}{procedure to do the inference}

  \item{G}{genetic relationship matrix}

  \item{...}{parameters passed to \code{genrel.matrix}}
}
\value{
numeric
}
\description{
Genome Scan With a Single QTL Model
}
\examples{
cross <- sim.cross.geno(250, nmar=10)
cross$pheno <- sim.cross.pheno(0.5, cross)
heritability(cross)
}
\keyword{manip}

