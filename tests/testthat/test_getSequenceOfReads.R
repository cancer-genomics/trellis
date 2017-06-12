context("Extract sequence of a read in a Rearrangement object")

test_that("getSequenceOfReads", {
  extdata <- system.file("extdata", package="svbams")
  id <- "cgov44t_revised.bam"
  bam.file <- file.path(extdata, id)
  data(rlist, package="svrearrange")
  improper_rp <- improper(rlist[[1]])

  with_seq <- DataFrame(getSequenceOfReads(rlist, bam.file, MAX=5L))
  expect_is(with_seq, "DataFrame")
  expect_identical(nchar(with_seq$seq[1]), 100L)
})
