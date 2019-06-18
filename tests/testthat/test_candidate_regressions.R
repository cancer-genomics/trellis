context("Regression tests for candidate rearrangements")

test_that("findCandidates", {
  library(svfilters.hg19)
  library(svbams)
  pdat <- cgov44t_preprocess()
  rparam <- RearrangementParams()
  rlist <- findCandidates2(pdat, rparam)
  if(FALSE){
    saveRDS(rlist, file="findCandidates.fe9b1f6.rds")
  }
  path <- system.file("extdata", package="svbams")
  rlist.fe9b1f6 <- readRDS(file.path(path, "findCandidates.fe9b1f6.rds"))
  expect_equivalent(rlist.fe9b1f6, rlist)
})

test_that("seqJunctionsInferredByPairedTags", {
  library(svfilters.hg19)
  pdat <- cgov44t_preprocess()
  rparam <- RearrangementParams()
  candidates <- seqJunctionsInferredByPairedTags2(pdat, rparam)
  if(FALSE){
    saveRDS(candidates, file="seqJunctionsInferredByPairedTags.fe9b1f6.rds")
  }
  path <- system.file("extdata", package="svbams")
  f <- file.path(path, "seqJunctionsInferredByPairedTags.fe9b1f6.rds")
  cand.fe9b1f6 <- readRDS(f)
  expect_equivalent(candidates, cand.fe9b1f6)
})

test_that("RearrangementList", {
  path <- system.file("extdata", package="svbams")
  f <- file.path(path, "seqJunctionsInferredByPairedTags.fe9b1f6.rds")
  candidates <- readRDS(f)
  cand.list <- RearrangementList(candidates)
  if(FALSE){
    saveRDS(cand.list, file="RearrangementList.fe9b1f6.rds")
  }
  path <- system.file("extdata", package="svbams")
  expected.fe9b1f6 <- readRDS(file.path(path, "RearrangementList.fe9b1f6.rds"))
  expect_equivalent(cand.list, expected.fe9b1f6)
})

test_that("type_each", {
  path <- system.file("extdata", package="svbams")
  cand.list <- readRDS(file.path(path, "RearrangementList.fe9b1f6.rds"))
  if(FALSE){
    for(j in seq_along(cand.list)){
      irp <- improper(cand.list[[j]])
      irp@first <- updateObject(irp@first)
      irp@last <- updateObject(irp@last)
      improper(cand.list[[j]]) <- irp
    }
  }
  result <- type_each(cand.list)
  if(FALSE){
    saveRDS(result, file="type_each.fe9b1f6.rds")
  }
  path <- system.file("extdata", package="svbams")
  result.fe9b1f6 <- readRDS(file.path(path, "type_each.fe9b1f6.rds"))
  expect_equivalent(result, result.fe9b1f6)
})

test_that("modalRearrangement", {
  path <- system.file("extdata", package="svbams")
  candidate.list <- readRDS(file.path(path, "type_each.fe9b1f6.rds"))
  colData(candidate.list)$modal_rearrangement <- modalRearrangement(candidate.list)
  colData(candidate.list)$percent_rearrangement <- percentRearrangement(candidate.list)
  path <- system.file("extdata", package="svbams")
  rlist.fe9b1f6 <- readRDS(file.path(path, "findCandidates.fe9b1f6.rds"))
  expect_identical(candidate.list, rlist.fe9b1f6)
})
