##context("Write improper alignments")
##test_that("writeImproperAlignments", {
####   library(svcontrols)
####   data(cgov2t_hg19)
####   writeImproperAlignments(cgov2t_hg19)
##
##  
##  ##  library(svovarian)
##  ##  dirs <- projectOvarian()
##  ##  bv <- readRDS(file.path(dirs[1], "bviews_hg19.rds"))
##  ##  test_views <- AlignmentViews2(bv, dirs)
##  ##  expect_is(test_views, "AlignmentViews2")
##  ##  expect_is(test_views[, 1], "AlignmentViews2")
##  ##  writeImproperAlignments2(test_views, mapq_thr=30)
##  ##  expect_true(all(file.exists(improperPaths(test_views))))
##  ##  ## this just reads from disk
##  ##  timing <- system.time(writeImproperAlignments2(test_views))
##  ##  expect_less_than(timing[["sys.self"]], 3)
##  ##
##  ##  gps <- readRDS(file.path(dirs["improper"], rdsId(test_views)[1]))
##  ##  expect_is(gps, "GAlignmentPairs")
##})

