test_that("organizeReads", {
  data(rlist, package="trellis")
  df.list <- organizeReads(rlist)
  expect_is(df.list, "list")
  expect_is(df.list[[1]], "data.frame")
  nms <- sapply(df.list, function(x) x$rid[1])
  expect_identical(nms, names(rlist))
})
