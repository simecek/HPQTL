\name{add.qtl}
\alias{add.qtl}
\title{Add a QTL to Simulated Phenotype}
\usage{
add.qtl(cross, qtl.pos, qtl.size, pheno.cols = 1)
}
\arguments{
  \item{cross}{"\code{cross}" object}

  \item{pheno.cols}{selection of phenotype's column(s)}

  \item{qtl.pos}{position of QTL(s), either number of
  markers or relative position in genome (between 0 and 1)}

  \item{qtl.size}{size of QTL(s),}
}
\value{
"\code{cross}" object
}
\description{
Add a QTL to Simulated Phenotype
}
\examples{
cross <- sim.cross.geno(250, nmar=10)
cross$pheno <- sim.cross.pheno(0.01, cross)
plot(scan1(cross))
cross <- add.qtl(cross, 0.75, 1)
plot(scan1(cross))
}
\keyword{manip}

