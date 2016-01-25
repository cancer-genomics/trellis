create_deletions_data <- function(){
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
  del_file <- file.path(dp["1deletions"], rdsId(bviews[, id]))
  if(file.exists(del_file)){
    file2 <- gsub(".rds", "_copy.rds", del_file)
    file.copy(del_file, file2)
    unlink(del_file)
  }
  grl <- readRDS(file.path(dp["segment"], "grl_hg19.rds"))
  gr <- grl[[id]]
  gr <- gr[gr$seg.mean < log2(0.75)]
  grl_del <- GRangesList(list(CGOV21T=gr))
  sv_dels <- sv_deletion_exp(dirs=dp, grl=grl_del[id],
                             bviews=bviews[, id],
                             gr_filters=germline_filters)
  deletions <- sv_dels[[1]]
  save(deletions, file="/dcl01/scharpf/data/svpackages/svclasses/data/deletions.rda")
  expect_is(sv_dels[[1]], "StructuralVariant")
  expect_true(file.exists(del_file))
  unlink(del_file)
}
