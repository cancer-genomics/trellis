context("Rearrangement")

test_that("improper<-", {
  data(rear_list, package="trellis")
  r <- rear_list[[1]]
  orig <- improper(r)
  ## Make sure the rear_list object doesn't change over time
  expect_equivalent(length(orig), 26)
  sub <- orig[1:2]
  improper(r) <- sub
  expect_equivalent(improper(r), sub)
})
