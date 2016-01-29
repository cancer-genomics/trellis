context("Call deletions")
test_that("deletions", {
  library(svovarian)
  library(svfilters)
  library(svcnvs)
  id <- "CGOV21T"
  dp <- projectOvarian(rootname="OvarianData2")
  data(lymphoblast_filters_hg19)
  data(lowMappabilityBins_hg19)
  data(binAssemblyGaps_hg19)
  gfilters <- as(lymphoblast_filters_hg19, "list")
  gfilters$map <- lowMappabilityBins_hg19
  gfilters$gc <- binAssemblyGaps_hg19
  germline_filters <- reduce(unlist(GRangesList(lapply(gfilters, granges))))
  bviews <- readRDS(file.path(dp[1], "bviews_hg19.rds"))
  grl <- readRDS(file.path(dp["segment"], "grl_hg19.rds"))
  gr <- grl[[id]]
  gr <- gr[gr$seg.mean < log2(0.75)]
  grl_del <- GRangesList(list(CGOV21T=gr))
  sv_dels <- sv_deletion_exp(dirs=dp, grl=grl_del[id],
                             bviews=bviews[, id],
                             gr_filters=germline_filters)
  expect_is(sv_dels[[1]], "StructuralVariant")
})

