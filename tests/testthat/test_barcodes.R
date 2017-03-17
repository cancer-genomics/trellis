context("Reading tags with barcodes")

test_that("readGAlignmentPairs", {
  library(GenomicAlignments)
  library(svbams)
  path <- system.file("extdata", package="svbams")
  targeted.bam <- file.path(path,
                            "PGDX5881P_PS_Seq_novo.sorted_nofa.bam")
  expect_true(file.exists(targeted.bam))
  gr <- GRanges("chr8", IRanges(5000, 100000))
  seqinfo(gr) <- Seqinfo(seqnames="chr8", seqlengths=146364022)
  flags <- scanBamFlag(isUnmappedQuery=FALSE,
                       hasUnmappedMate=FALSE,
                       isSecondaryAlignment=FALSE)
  param <- ScanBamParam(flag=flags,
                        ##which=gr,
                        mapqFilter=30,
                        tag="BC")
  galp <- readGAlignmentPairs(targeted.bam, param=param)
})
