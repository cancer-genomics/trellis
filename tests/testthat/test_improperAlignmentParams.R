context("improper alignment parameters")

test_that("improper alignment parameters", {
  library(Rsamtools)
  params <- improperAlignmentParams()
  expect_true(bamFlag(params)[["isPaired"]])
  expect_true(!bamFlag(params)[["isProperPair"]])
  expect_true(!bamFlag(params)[["isUnmappedQuery"]])
  expect_true(!bamFlag(params)[["hasUnmappedMate"]])
  expect_true(!bamFlag(params)[["isDuplicate"]])
  expect_true(is.na(bamFlag(params)[["isNotPassingQualityControls"]]))
  flags <- improperAlignmentFlags()
  params2 <- improperAlignmentParams(flag=flags)
  expect_identical(params, params2)

  params3 <- improperAlignmentParams(flag=flags, mapqFilter=30)
  expect_identical(bamMapqFilter(params3), 30L)
  params3 <- improperAlignmentParams(flag=flags, mapqFilter=0)

  params <- properAlignmentParams()
  expect_true(bamFlag(params)[["isPaired"]])
  expect_true(bamFlag(params)[["isProperPair"]])
  expect_true(!bamFlag(params)[["isUnmappedQuery"]])
  expect_true(!bamFlag(params)[["hasUnmappedMate"]])
  expect_true(!bamFlag(params)[["isDuplicate"]])
  expect_true(is.na(bamFlag(params)[["isNotPassingQualityControls"]]))
  flags <- properAlignmentFlags()
  pparams2 <- improperAlignmentParams(flag=flags, mapqFilter=30)
  expect_identical(bamMapqFilter(pparams2), 30L)

  pparams2 <- improperAlignmentParams(flag=flags, mapqFilter=0)
  expect_identical(bamMapqFilter(pparams2), 30L)
})

.test_that <- function(expr, ...) NULL

test_that("read_pairs_from_bam", {
  library(Rsamtools)
  library(svbams)
  library(svalignments)
  path <- system.file("extdata", package="svbams")

  library(TxDb.Hsapiens.UCSC.hg19.refGene)
  region <- GRanges("chr15", IRanges(63201003-2000, 63209243+2000))
  si <- seqinfo(TxDb.Hsapiens.UCSC.hg19.refGene)
  seqinfo(region) <- si["chr15", ]

  bampath <- list.files(path, pattern="cgov44t.bam$", full.names=TRUE)
  bview <- BamViews(bamPaths=bampath)

  iparams <- improperAlignmentParams(which=region, mapqFilter=30)
  pparams <- properAlignmentParams(which=region, mapqFilter=30)
  irp <- getImproperAlignmentPairs(bview,
                                   iparams)
  expect_identical(length(irp), 57L)
  g.irp <- ga2gr(irp, is.improper=TRUE)
  expect_identical(length(g.irp), 57L*2L)
  expect_true(all(g.irp$is.improper))

  prp <- getProperAlignmentPairs(bview,
                                 pparams)
  g.prp <- ga2gr(prp, is.improper=FALSE)
  mcol.vars <- colnames(mcols(g.prp))
  expect_identical(mcol.vars, c("read", "is.improper", "tagid"))
  expect_identical(length(g.prp), length(prp) * 2L)
  expect_true(!any(g.prp$is.improper))
  L <- length(g.prp)
  g.prp2 <- thinProperPairs(g.prp, 10)
  expect_equal(length(g.prp2), length(g.prp)/10, tolerance=5)
  gr <- c(g.irp, g.prp)
  ## Because reads are paired, half the length should be the number of
  ## duplicated tags
  expect_identical(length(gr)/2, sum(duplicated(gr$tagid))*1)
  gr <- sortByRead1(gr)

  ##
  ## ga2gr should return a length-0 GRanges object when there are no
  ## proper/improper alignment pairs
  ##
  prp2 <- prp[rep(FALSE, length(prp))]
  lengthzero.gr <- ga2gr(prp2, is.improper=FALSE)
  expect_is(lengthzero.gr, "GRanges")
  expect_identical(length(lengthzero.gr), 0L)
})






