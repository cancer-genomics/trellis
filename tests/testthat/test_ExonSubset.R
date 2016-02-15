context("ExonSubset")
test_that("ExonSubset", {
  es <- ExonSubset()
  expect_true(validObject(es))

})
