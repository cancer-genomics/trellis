context("writeToFasta")

test_that("writeToFasta", {
  path <- system.file("extdata", package = "trellis")
  tags <- readRDS(file.path(path, "tags.34bd277.rds"))
  fa.file <- tempfile()
  writeToFasta(tags, fa.file)
  if(FALSE){
    writeToFasta(tags, "fasta_34bd277.txt")
  }
  fa <- read.delim(fa.file)
  expected <- read.delim("fasta_34bd277.txt")
  expect_identical(fa, expected)
})

test_that("umapped_read", {
  path <- system.file("extdata", package = "trellis")
  rlist <- readRDS(file.path(path, "filterRearrangementList.37be188.rds"))
  query <- uncouple(linkedBins(rlist))
  extdata <- system.file("extdata", package="svbams")
  bam.file <- file.path(extdata, "cgov44t_revised.bam")
  unmapped <- unmapped_read(bam.file, query, yield_size=200000)
  if(FALSE){
    saveRDS(unmapped, file="unmapped.34bd277.rds")
  }
  expected <- readRDS(file.path(path, "unmapped.34bd277.rds"))
  expect_identical(expected, unmapped)
})
