% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/sim.cross.pheno.r
\name{sim.cross.pheno}
\alias{sim.cross.pheno}
\title{Generate Random Phenotype with Given Heritability}
\usage{
sim.cross.pheno(h2, geno, background = "GRM", ...)
}
\arguments{
\item{h2}{heritability}

\item{geno}{"\code{genotype.probs}" object}

\item{background}{method to generate genetic background, see Details}

\item{...}{parameters passed to \code{genrel.matrix}}
}
\value{
matrix with phenotype data
}
\description{
Generate Random Phenotype with Given Heritability
}
\details{
Currently two polygenic backgrounds are supported:
If \code{background = "GRM"} then the genetic effect is generated from a multivariate
normal distribution with genetic similarity matrix as variance matrix. If
\code{background = "all-snps"} then genetic effect is a sum of small gaussian distributed
effect at every SNP location.
}
\examples{
cross <- sim.cross.geno(250, nmar=10)
cross$pheno <- sim.cross.pheno(0.5, cross)
cross$pheno <- sim.cross.pheno(0.5, cross, method = "kinship")
cross$pheno <- sim.cross.pheno(0.5, cross, background = "all-snps")
}
\keyword{manip}

