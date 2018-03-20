context("properReadPairs")

test_that("properReadPairs",{
  library(svbams)
  bamdir <- system.file("extdata", package="svbams", mustWork=TRUE)
  bam.file <- file.path(bamdir, "cgov10t.bam")
  ##svfile <- file.path(bamdir, "cgov10t_deletions.rds")
  bins <- readRDS(file.path(bamdir, "cgov10t_bins1kb.rds"))
  proper <- properReadPairs(bam.file, bins)

  ## does not do anything if all bins are less than 20kb
  bins2 <- thinReadPairQuery(bins)
  expect_equivalent(bins2, bins)

  segments <- variant(readRDS(file.path(bamdir, "cgov10t_deletions.rds")))
  ## for larger segments, focus on the border 
  segs <- thinReadPairQuery(segments)
  expect_identical(length(segs), 10L)
})
