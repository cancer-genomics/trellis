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

test_that("rbind", {
  library(S4Vectors)
  df1 <- DataFrame(a=letters)
  df2 <- DataFrame(a=letters)
  result <- rbind(df1, df2)
  expect_is(result, "DataFrame")
})

test_that("c", {
  tmp <- c(RearrangementList(), RearrangementList())
  expect_is(tmp, "RearrangementList")
  data(rear_list)
  tmp <- c(rear_list[1:2], rear_list[3:4])
  expect_is(tmp, "RearrangementList")
  expect_equivalent(tmp, rear_list[1:4])
  ##
  ## NOTE: for some reason do.call(rbind, list_of_DataFrames) does not work
  ##  - had to coerce to data.frame and then back to DataFrame
  ##  - this causes the colData slot not to be identical
})

test_that("splitReads<-", {
  data(rear_list)
  nms <- names(rear_list)
  g <- replicate(4, GRanges())
  grl <- GRangesList(g)
  names(grl) <- head(nms, 4)
  splitReads(rear_list) <- grl

  ##
  ## accessor
  ##
  sr <- head(splitReads(rear_list), 4)
  expect_equivalent(sr, grl)
})


