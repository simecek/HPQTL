\name{gensim.matrix}
\alias{gensim.matrix}
\title{Genetic Similarity Matrix}
\usage{
gensim.matrix(geno, method = c("default", "allele-2f-additive",
  "allele-multif-additive", "allele-multif-cosine"), procedure = c("LMM",
  "LMM-L1O", "LM"), subjects = seq(geno$subjects),
  markers = seq(NROW(geno$markers)))
}
\arguments{
  \item{geno}{\code{genotype.probs} object or \code{cross}
  object}

  \item{method}{how to calculate similarity matrix
  ('allele-2f-additive', 'allele-multif-additive')}
}
\value{
square matrix (\code{nind} rows, \code{nind} columns)
}
\description{
Genetic Similarity Matrix
}
\examples{
cross <- sim.cross.geno(250, nmar=10)
geno <- extract.geno(cross)
gensim.matrix(geno)
gensim.matrix(geno, method = "allele-multif-additive")
}
\keyword{manip}

