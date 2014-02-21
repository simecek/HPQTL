\name{add.qtl}
\alias{add.qtl}
\title{Add a QTL to Simulated Phenotype}
\usage{
add.qtl(cross, qtl.pos, qtl.size, pheno.cols = 1)
}
\arguments{
  \item{cross}{"\code{cross}" object}

  \item{pheno.cols}{selection of phenotype's column(s)}

  \item{qtl.pos}{position of QTL(s), qither number of
  markets or relative position in genome (between 0 and 1)}

  \item{qtl.size}{size of QTL(s),}
}
\value{
numeric
}
\description{
Add a QTL to Simulated Phenotype
}
\examples{
cross <- sim.cross.geno(250, nmar=10)
cross <- calc.genoprob(cross)
plot(scanone(cross))
cross.with.qtls <- add.qtl(cross, c(0.15,0.35), c(1,1))
plot(scanone(cross.with.qtls))
}
\keyword{manip}

