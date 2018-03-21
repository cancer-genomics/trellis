context("complex rearrangement")

test_that("complex", {
  extdata <- system.file("extdata", package="svrearrange")
  r <- readRDS(file.path(extdata, "cgov1t_complex_rearrangement.rds"))
  s <- type(r)
  expect_identical(s, "++,--")
  expect_true(isComplex(r))
})
