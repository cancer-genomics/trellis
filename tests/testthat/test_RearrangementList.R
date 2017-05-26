context("RearrangementList")

test_that("Rearrangement", {
  Rearrangement2 <- function(linkedBins=GRanges(),
                             improper=.GAlignmentPairs(),
                             partitioning=integer(),
                             link1id=character(),
                             link2id=character(),
                             tags=GRanges(),
                             tag_map_to_linked_bin=character(),
                             modal_rearrangement=character(),
                             percent_rearrangement=numeric(),
                             fraction_linking_tags=numeric(),
                             split_reads=GRanges()){
    new("Rearrangement",
        linkedBins=linkedBins,
        improper=improper,
        partitioning=partitioning,
        link1id=link1id,
        link2id=link2id,
        tags=tags,
        tag_map_to_linked_bin=tag_map_to_linked_bin,
        modal_rearrangement=modal_rearrangement,
        percent_rearrangement=percent_rearrangement,
        fraction_linking_tags=fraction_linking_tags,
        split_reads=GRanges())
  }
  r <- Rearrangement2()
  expect_true(validObject(r))
  expect_is(splitReads(r), "GRanges")
  splitReads(r) <- GRanges()
  expect_is(splitReads(r), "GRanges")
})

test_that("RearrangementList", {
  ## constructor
  rl <- RearrangementList()
  expect_true(validObject(rl))
  rl$something <- character()
  expect_identical(rl$something, character())
  expect_identical(nrow(colData(rl)), 0L)
  expect_identical(ncol(colData(rl)), 1L)
})
