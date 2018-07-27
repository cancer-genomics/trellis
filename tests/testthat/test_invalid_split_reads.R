context("Spit read checks")

test_that("is_valid_splits", {
  extdata <- system.file("extdata", package="svbams")
  rlist <- readRDS(file.path(extdata, "rlist_endometrioid_project.rds"))
  if(FALSE){
    for(i in seq_along(rlist)){
      r <- rlist[[i]]
      irp <- improper(r)
      first(irp) <- updateObject(irp@first)
      last(irp) <- updateObject(irp@last)
      improper(r) <- irp
      rlist[[i]] <- r
    }
    saveRDS(rlist, file.path(extdata, "rlist_endometrioid_project.rds"))
  }
  expect_true(all(is_valid_splits(rlist, maxgap=200)))
  expect_false(all(is_valid_splits(rlist, maxgap=50)))
  n.sr1 <- elementNROWS(splitReads(rlist))
  rlist2 <-  fiveTo3List(rlist, "hg18")
  ## two orientations for each rearrangement
  expect_equal(length(rlist2), length(rlist)*2)
  n.sr2 <- elementNROWS(splitReads(rlist2))
  ## note that some of the split reads have been removed
  ## when constructing rearrangements with a specific orientation
  ## - this happens because one of the rearrangements contains
  ## split reads with both ++ and -- strands
  expect_true(any(n.sr2[names(n.sr1)] < n.sr1))
  ## now some of the rearrangements are no longer supported
  ## by split reads even with the more lenient gap
  is_valid <- is_valid_splits(rlist2, maxgap=200)
  rlist2 <- rlist2[ is_valid ]

  ## using a more stringent maxgap may reduce the
  ## need for filtering after fiveTo3List
  rlist <- readRDS(file.path(extdata, "rlist_endometrioid_project.rds"))
  is_valid <- is_valid_splits(rlist, maxgap=50)
  rlist <- rlist[ is_valid ]
  rlist2 <-  fiveTo3List(rlist, "hg18")
  is_valid <- is_valid_splits(rlist2)
  expect_true(all(is_valid))
  ##
  ## Test for exception when no split reads present
  ##
  splitReads(rlist2[[1]]) <- splitReads(rlist2[[1]])[0]
  is_valid <- is_valid_splits(rlist2[1])
  ## not sure why we want to allow this to be valid. Maybe we should require at least 1 split read here
  expect_true(is_valid)
  rlist2 <- fiveTo3List(rlist2[1], build="hg18")
  expect_error(rearDataFrameList(rlist2))
})
