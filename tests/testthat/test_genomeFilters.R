context("Germline filters")
test_that("germlineFilters", {
  gfs <- listGenomeFilters()
  expect_is(gfs, "list")
})
