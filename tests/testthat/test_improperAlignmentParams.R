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

  irp <- getImproperAlignmentPairs(bview,
                                   iparams,
                                   mapq_thr=30,
                                   use.mcols=TRUE)
  prp <- getProperAlignmentPairs(bview,
                                 pparams,
                                 mapq_thr=30,
                                 use.mcols=TRUE)
  readPairs <- function(bview, mapq_thr=30){
    irp <- getImproperAlignmentPairs(bview,
                                     iparams,
                                     mapq_thr=30,
                                     use.mcols=TRUE)
    prp <- getProperAlignmentPairs(bview,
                                   pparams,
                                   mapq_thr=30,
                                   use.mcols=TRUE)
    ##mcols(prp)$is.improper <- FALSE
    ##mcols(irp)$is.improper <- TRUE
  }

  r1.irp <- granges(first(irp))
  r2.irp <- granges(last(irp))
  names(r1.irp) <- names(r2.irp) <- NULL
  r1.irp$read <- "R1"
  r2.irp$read <- "R2"
  r1.irp$is.improper <- r2.irp$is.improper <- TRUE
  r1.prp <- granges(first(prp))
  r2.prp <- granges(last(prp))
  names(r1.prp) <- names(r2.prp) <- NULL
  r1.prp$read <- "R1"
  r2.prp$read <- "R2"
  r1.prp$is.improper <- r2.prp$is.improper <- FALSE
  r1 <- c(r1.prp, r1.irp)
  r2 <- c(r2.prp, r2.irp)
  r1$pair.id <- seq_len(length(r1))
  r2$pair.id <- seq_len(length(r2))
  g <- c(r1, r2)
})




