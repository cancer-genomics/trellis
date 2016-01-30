context("RearrangementList")
test_that("RearrangementList", {
  ## constructor
  rl <- RearrangementList()
  expect_true(validObject(rl))
  rl$something <- character()
  expect_identical(rl$something, character())
  expect_identical(nrow(colData(rl)), 0L)
  expect_identical(ncol(colData(rl)), 1L)
})
