context("seqinfo")
test_that("seqinfo", {
  pv <- PreprocessViews2()
  expect_identical(seqlevels(pv), character())
  expect_is(seqinfo(pv), "Seqinfo")
  seqinfo(pv) <- Seqinfo()
})
