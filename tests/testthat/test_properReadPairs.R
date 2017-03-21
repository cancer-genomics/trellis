context("properReadPairs")

test_that("properReadPairs",{
  library(svbams)
  svfile <- file.path(bamdir, "cgov10t_deletions.rds")
  bamdir <- system.file("extdata", package="svbams")
  bins <- readRDS(file.path(bamdir, "cgov10t_bins1kb.rds"))
  bview <- BamViews(bamPaths=bamfile, bamRanges=bins)
})
