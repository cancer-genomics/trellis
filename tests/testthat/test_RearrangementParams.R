context("Rearrangement Parameters")
test_that("rearrangementParams", {
  rp <- RearrangementParams()
  expect_true(validObject(rp))

  ## accessors
  ## minimum separation of read pairs to be considered 'improper'
  expect_identical(rpSeparation(rp), 10e3)

  ## minimum number of single tags required to constitute a 'cluster'
  expect_identical(minNumberTagsPerCluster(rp), 5)

  ## if we take the union of the tags that form a cluster, their width
  ## must be at least minClustersize() basepairs
  expect_identical(minClusterSize(rp), 115L)

  ## the size of the cluster can not be larger than maxClusterSize()
  expect_identical(maxClusterSize(rp), 5000L)

  ## tags with a gap of less than minGapWidth() are considered part of
  ## the same cluster
  expect_identical(minGapWidth(rp), 1000L)

  ## the fraction of linking tags that support the modal rearrangement
  expect_identical(percentModalType(rp), 0.9)

  ## at least percentModalType() of the linking tags must support the
  ## modal rearrangement
  expect_identical(percentLinking(rp), 0.8)


  rp2 <- RearrangementParams(min_cluster_size=115L,
                             max_cluster_size=5e3L,
                             min.gapwidth=1e3L)
  expect_identical(rp2, RearrangementParams())
})
