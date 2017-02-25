test_that("DeletionParam", {
  dp <- DeletionParam()
  expect_true(validObject(dp))

  dp <- DeletionParam(max_tag_density=2)
  homozygousThr(dp)
})
