\name{scan1}
\alias{scan1}
\title{Genome Scan With a Single QTL Model}
\usage{
scan1(pheno, geno, covar = NULL, procedure = c("LM", "LMM", "LMM-L1O"),
  G = NULL, ...)
}
\arguments{
  \item{cross}{"\code{cross}" object}

  \item{pheno.cols}{selection of phenotype's column(s)}

  \item{covar}{(additive) covariates}

  \item{procedure}{procedure to do the inference, see
  Details}

  \item{V}{genetic similarity matrix or a list of genetic
  similarity matrices or \code{NULL}}

  \item{...}{parameters passed to \code{genrel.matrix}}
}
\value{
\code{scanone} object
}
\description{
Genome Scan With a Single QTL Model
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
cross$pheno <- sim.cross.pheno(c(0.01, 0.5, 0.75), cross)
plot(scan1(cross, pheno.cols=1:3), lodcol=1:3)
plot(scan1(cross, pheno.cols=1:3, procedure="scanOne"), lodcol=1:3)

# slow - do not run
# plot(scan1(cross, pheno.cols=1:3, procedure="scanOne-per-chr"), lodcol=1:3)
}
\keyword{manip}

