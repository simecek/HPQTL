\name{scan1.threshold}
\alias{scan1.threshold}
\title{Significance Threshold of Genome Scan}
\usage{
scan1.threshold(cross, pheno.cols = 1, procedure = "scanone",
  n.perm = 100, alpha = 0.05, ...)
}
\arguments{
  \item{cross}{"\code{cross}" object}

  \item{pheno.cols}{selection of phenotype's column(s)}

  \item{procedure}{procedure to do the inference}

  \item{n.perm}{number of permutations}

  \item{alpha}{level of significance}

  \item{...}{parameters passed to \code{genrel.matrix}}
}
\value{
numeric
}
\description{
Significance Threshold of Genome Scan
}
\examples{
cross <- sim.cross.geno(250, nmar=10)
cross$pheno <- sim.cross.pheno(0.5, cross)
heritability(cross)
}
\keyword{manip}

