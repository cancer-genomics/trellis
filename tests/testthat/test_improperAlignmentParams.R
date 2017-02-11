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
})

.test_that <- function(expr, ...) NULL

.test_that("read_pairs_from_bam", {
  library(Rsamtools)
  library(svbams)
  library(svalignments)
  path <- system.file("extdata", package="svbams")

  library(TxDb.Hsapiens.UCSC.hg19.refGene)
  region <- GRanges("chr15", IRanges(63201003, 63209243))
  si <- seqinfo(TxDb.Hsapiens.UCSC.hg19.refGene)
  seqinfo(region) <- si["chr15", ]

  bampath <- list.files(path, pattern="cgov44t.bam$", full.names=TRUE)
  bview <- BamViews(bamPaths=bampath)

  iparams <- improperAlignmentParams()
  pparams <- properAlignmentParams()
  irp <- getImproperAlignmentPairs(bview,
                                   iparams,
                                   mapq_thr=30,
                                   use.mcols=TRUE)
  expect_identical(length(irp), 445L)
  g.irp <- ga2gr(irp, is.improper=TRUE)
  expect_identical(length(g.irp), 445L*2L)
  expect_true(all(g.irp$is.improper))

  prp <- getProperAlignmentPairs(bview,
                                 pparams,
                                 mapq_thr=30,
                                 use.mcols=TRUE)
  g.prp <- ga2gr(prp, is.improper=FALSE)
  expect_identical(length(g.prp), length(prp) * 2L)
  expect_true(!any(g.prp$is.improper))
})




