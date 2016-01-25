context("DataPaths")
test_that("DataPaths", {
  dp <- DataPaths()

  ## expect the first directory to be data
  expect_identical(basename(svclasses:::dataDir(dp)), "data")
  
  preprocess_dirs <- listDir("preprocess", dp)
  expect_identical(length(preprocess_dirs[["folders"]]), 5L)

  segment_dirs <- listDir("segment", dp)
  expect_identical(length(segment_dirs[["folders"]]), 4L)

  projects <- topics(dp)
})
