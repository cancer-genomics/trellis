test_that("axis_limits", {
  extdir <- system.file("extdata", package="trellis")
  rear <- readRDS(file.path(extdir, "chr6-12_rear.rds"))
  tmp <- .rear_DataFrame(rear)
  expect_true(tmp$junction_3p[1] == 24259912)

  df <- rearDataFrame(rear)
  limits <- axis_limits(df, 400)
  sizes <- sapply(limits, diff) %>%
    abs
  expect_true(all(sizes < 500))
})
