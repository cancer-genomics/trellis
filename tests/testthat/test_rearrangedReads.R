test_that("overlapsBoth", {
  data(rearrangement_list)
  data(blat_unmapped)
  data(rearrangement_list)
  rearranged_reads <- rearrangedReads(rearrangement_list, blat_unmapped, 500)
  expect_is(rearranged_reads, "GRangesList")
  ## only one rearrangement supported
  expect_identical(length(rearranged_reads), 1L)
  expect_true(length(rearrangement_list[names(rearranged_reads)])==1)
  ## number split reads
  n_splitreads <- length(split(rearranged_reads[[1]], rearranged_reads[[1]]$qname))
  expect_identical(n_splitreads, 32L)
  
  ## this should exit gracefully with no rearranged reads identified
  rearranged_reads <- rearrangedReads(rearrangement_list, blat_unmapped[1:10, ], 500)
  expect_identical(length(rearranged_reads), 0L)
  expect_is(rearranged_reads, "GRangesList")
})
