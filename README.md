HPQTL
=====

HPQTL is R package for quantitative trait locus (QTL) analysis in case of highly polygenic genetic backgound.

To calculate LOD scores, the user need to specify the following
* phenotype (variable of interest + covariates)
* genotype (3-dim array of probabilities: subjects x calls x markers)
* genetic similarity matrix (can be calculated from genotype)

Three methods have been implemented
* linear model (LM)
* linear mixed model (LMM)
* linear mixed model with specific genetic similarity matrix for each chromosome (LMM-L1O)

## Example

```S
library(HPQTL)
data(fake.f2, package="qtl")

# calculate and extract genotype probabilities
fake.f2 <- qtl::calc.genoprob(fake.f2)
geno <- extract.geno(fake.f2)

# calculate genetic similarity matrix
G <- gensim.matrix(geno)

# mapping with linear model
qtl::scanone(fake.f2, method = "hk")
(fit.lm <- scan1(fake.f2$pheno[,1], geno=geno))

# mapping with linear mixed model
fit.lmm <- scan1(fake.f2$pheno[,1], geno=geno, procedure = "LMM", G=G)

# mapping with linear mixed model - leave the scanned chromosome out
fit.lmm1 <- scan1(fake.f2$pheno[,1], geno=geno, procedure = "LMM-L1O")

# LOD plots
plot(fit.lm, col="black")
plot(fit.lmm, add=TRUE, col="red")
plot(fit.lmm1, add=TRUE, col="blue")
legend("topleft", c("LM", "LMM", "LMM-L1O"), lty=1, col=c("black", "red", "blue"))

```

## Installation

To install the package directly from Github, use `devtools::install_github`:

```S
    library(devtools)
    install_github("simecek/HPQTL")
```

## See also

* [QTLRel](https://github.com/pcarbo/QTLRel)
* [emma](http://mouse.cs.ucla.edu/emma/)
* [qtl](https://github.com/kbroman/qtl)
* [qtlcharts](https://github.com/kbroman/qtlcharts)

