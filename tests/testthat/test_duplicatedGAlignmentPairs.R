test_that("duplicatedGAlignmentPairs", {
  library(Rsamtools)
  library(svbams)
  library(GenomicAlignments)
  data(bins)
  ### Read in bam file and great GAlignmentPairs object
  bam_path <- system.file("extdata", package="svbams")
  bamfile <- file.path(bam_path, "cgov10t.bam")
  data(bins, package='svbams')
  param1 <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE,
                                            isSecondaryAlignment = FALSE),
                         which=bins)
  param1b <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE,
                                             isSecondaryAlignment = FALSE))
  galp1 <- readGAlignmentPairs(bamfile, param=param1)
  galp1 <- galp1[do.call(order, as.data.frame(galp1))]
  galp1 <- galp1[!duplicatedGAlignmentPairs(galp1)]
  galp1b <- readGAlignmentPairs(bamfile, param=param1b)
  galp1b <- subsetByOverlaps(galp1b, bins, ignore.strand=TRUE)
  galp1b <- galp1b[do.call(order, as.data.frame(galp1b))]
  answer <- TRUE 
  result <- identical(galp1, galp1b)
  expect_identical(result, answer) 
})
