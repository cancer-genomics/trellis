context("Create table of fusions")

test_that("fusion_table", {
  library(GenomicRanges)
  library(svfilters.hg19)
  extdata <- system.file("extdata", package="svbams")
  fusions <- readRDS(file.path(extdata, "valid_fusions.rds"))
  coding_jxns <- readRDS(file.path(extdata, "coding_jxns.rds"))
  ##
  ## hgnc symbols
  ##
  df <- fusionTable2(fusions)
  expect_is(df, "DataFrame")
})
