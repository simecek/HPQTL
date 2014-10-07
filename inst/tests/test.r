
context('LOD scores equal for HPQTL/scan1/LM and qtl/scanone/hk')

test_that("fake.f2",{

  data(fake.f2, package="qtl")
  fake.f2 <- calc.genoprob(fake.f2)
  
  f2.qtl <- scanone(fake.f2, method="hk")
  f2.hpqtl.lm <- scan1(fake.f2, procedure="LM")
  
  # expect qtl::scanone and scan1(..., method="LM") results equal
  expect_equal(f2.qtl$lod[f2.qtl$chr!="X"], f2.hpqtl.lm$lod[f2.hpqtl.lm$chr!="X"])
})

test_that("fake.bc",{

  data(fake.bc, package="qtl")
  fake.bc <- calc.genoprob(fake.bc)
  
  bc.qtl <- scanone(fake.bc, method="hk")
  bc.hpqtl.lm <- scan1(fake.bc, procedure="LM")
  
  # expect qtl::scanone and scan1(..., method="LM") results equal
  expect_equal(bc.qtl$lod[bc.qtl$chr!="X"], bc.hpqtl.lm$lod[bc.hpqtl.lm$chr!="X"])
})

context('LOD scores equal for HPQTL/scan1/LMM and qtlRel/scanOne')

test_that("fake.f2",{

  data(fake.f2, package="qtl")
  fake.f2 <- calc.genoprob(fake.f2)
  f2.qtl <- scanone(fake.f2, method="hk")
  
  G <- gensim.matrix(fake.f2)
  f2.hpqtl.lmm <- scan1(fake.f2, procedure="LMM", package="QTLRel")
  
  prDat <- list()
  prDat$pr <- extract.geno(fake.f2)$probs
  prDat$chr <- f2.qtl$chr[f2.qtl$chr!="X"]
  prDat$dist <- f2.qtl$pos[f2.qtl$chr!="X"]
  prDat$snp <- rownames(f2.qtl)[f2.qtl$chr!="X"]
  vc <- estVC(y=fake.f2$pheno[,1],v=list(AA=G, DD=NULL, HH=NULL, AD=NULL,MH=NULL,EE=diag(nrow(G))))
  f2.qtlrel <- scanOne(y=fake.f2$pheno[,1], x=rep(1,nrow(fake.f2$pheno)), gdat=NULL, prdat=prDat, vc=vc)$p / (2*log(10))

  # expect QTLRel::scanOne and scan1(..., method="LMM") results equal
  expect_equal(as.vector(f2.qtlrel), as.vector(f2.hpqtl.lmm$lod))
})

context('LOCO results != LMM results')

test_that("fake.f2",{
  
  data(fake.f2, package="qtl")
  fake.f2 <- calc.genoprob(fake.f2)

  f2.hpqtl.lmm <- scan1(fake.f2, procedure="LMM")
  f2.hpqtl.lmm_l1o <- scan1(fake.f2, procedure="LOCO")
  expect_false(all(f2.hpqtl.lmm$lod == f2.hpqtl.lmm_l1o$lod))
})

#context('Permutation thresholds')

#test_that("fake.f2",{
#  
#  data(fake.f2, package="qtl")
#  fake.f2 <- calc.genoprob(fake.f2)
#  
#  trhold.lm <- scan1.threshold(fake.f2, procedure="LM")
#  trhold.lmm <- scan1.threshold(fake.f2, procedure="LMM")
#  trhold.lmm_l1o <- scan1.threshold(fake.f2, procedure="LOCO")
#})
