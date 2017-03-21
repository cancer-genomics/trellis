context("Read improperly paired reads")

test_that("getImproperAlignmentPairs", {
  library(Rsamtools)
  library(svbams)
  bamdir <- system.file("extdata", package="svbams")
  bamfile <- file.path(bamdir, "cgov10t.bam")
  cgov10t.bed <- read.delim("../../../svbams/data-raw/cgov10t.bed", header=FALSE)
  chrom <- cgov10t.bed$V1
  g <- GRanges(chrom, IRanges(cgov10t.bed$V2, cgov10t.bed$V3))
  params <- improperAlignmentParams(which=g, mapqFilter=0)
  ##trace(getImproperAlignmentPairs, browser)
  ## 476 alignments with ambiguous pairing dumped
  irp <- getImproperAlignmentPairs(bamfile, param=params)
  sequence.names <- names(table(c(chromosome(first(irp)), chromosome(last(irp)))))
  expect_identical(sequence.names, unique(chrom))
})
