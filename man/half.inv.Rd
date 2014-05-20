\name{half.inv}
\alias{half.inv}
\title{Half-inverse matrix operation}
\usage{
half.inv(W, symmetric = TRUE, inverse = TRUE)
}
\arguments{
  \item{W}{matrix}
}
\description{
Half-inverse matrix operation
}
\examples{
data(fake.f2, package="qtl")
fake.f2 <- calc.genoprob(fake.f2)

geno <- extract.geno(fake.f2)
G <- gensim.matrix(geno)
A <- HPQTL:::half.inv(G)
}
\author{
copied from QTLRel
}
\keyword{manip}

