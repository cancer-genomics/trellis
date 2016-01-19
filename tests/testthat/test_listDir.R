context("DataPaths")
test_that("DataPaths", {
  dp <- DataPaths()
  preprocess_dirs <- listDir("preprocess", dp)
  expect_identical(length(preprocess_dirs[["folders"]]), 4L)

  segment_dirs <- listDir("segment", dp)
  expect_identical(length(segment_dirs[["folders"]]), 3L)  
})
