\name{half.inv}
\alias{half.inv}
\title{W^(-1/2) matrix operation}
\usage{
half.inv(W, symmetric = TRUE, inverse = TRUE)
}
\arguments{
  \item{W}{matrix}
}
\value{
matrix \code{A} such that \code{t(A)%*%A = W^{-1}}
}
\description{
W^(-1/2) matrix operation
}
\examples{
cross <- sim.cross.geno(250, nmar=10)
cross$pheno <- sim.cross.pheno(0.5, cross)
heritability(cross)
}
\author{
copied from QTLRel
}
\references{
Abney M, Ober C, McPeek MS (2002). "Quantitative trait
homozygosity and association mapping and empirical
genome-wide significance in large complex pedigrees:
Fasting serum insulin levels in the Hutterites." American
Journal of Human Genetics 70: 920-934.
}
\keyword{manip}

