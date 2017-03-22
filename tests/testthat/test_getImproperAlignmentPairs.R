context("Read improperly paired reads")

test_that("getImproperAlignmentPairs", {
  library(Rsamtools)
  library(svbams)
  library(GenomicAlignments)
  bamdir <- system.file("extdata", package="svbams")
  bamfile <- file.path(bamdir, "cgov10t.bam")
  params <- improperAlignmentParams(mapqFilter=0)
  irp <- getImproperAlignmentPairs(bamfile, param=params)
  sequence.names <- names(table(c(chromosome(first(irp)), chromosome(last(irp)))))
  expect_identical(sequence.names, c("chr3", "chr4", "chr9", "chrX"))
  ## add one regression test
  expect_identical(length(irp), 882L)
})
