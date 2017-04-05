test_that("Transcript", {
  cl.tx <- ClippedTranscripts()
  expect_true(validObject(cl.tx))
  nms <- c("5prime", "3prime")
  expect_identical(names(cl.tx), nms)
})
