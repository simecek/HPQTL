\name{scan1.threshold}
\alias{scan1.threshold}
\title{Significance Threshold of Genome Scan}
\usage{
scan1.threshold(pheno, geno, covar = NULL, procedure = c("LM", "LMM",
  "LMM-L1O"), G = NULL, n.perm = 100, alpha = 0.05, ...)
}
\arguments{
  \item{cross}{"\code{cross}" object}

  \item{pheno.cols}{selection of phenotype's column(s)}

  \item{procedure}{procedure to do the inference, see
  Details}

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
\details{
Currently, three procedures fully implemented:
\code{procedure = "qtl"} calls \code{scanone} function from
\code{qtl} package, \code{procedure = "QTLRel"} calls
\code{scanOne} function from \code{QTLRel} package,
\code{procedure = "QTLRel-pre-chr"} calls \code{scanOne}
for each chromosome separately and genetic relationship
matrix is estimated with this chromosome excluded.
}
\examples{
cross <- sim.cross.geno(250, nmar=10)
cross$pheno <- sim.cross.pheno(0.5, cross)
scan1.threshold(cross)

# slow - do not run
# scan1.threshold(cross, n.perm=1000, alpha=0.01)
# scan1.threshold(cross, procedure = "scanOne")
# scan1.threshold(cross, procedure = "scanOne-per-chr")
}
\keyword{manip}

