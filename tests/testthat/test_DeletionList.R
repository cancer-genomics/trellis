context("DeletionList")
test_that("DeletionList", {
  dl <- DeletionList()
  expect_true(validObject(dl))
})
