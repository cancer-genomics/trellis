context("AmpliconGraph")

test_that("AmpliconGraph", {
  ag <- AmpliconGraph()
  expect_true(validObject(ag))

  library(svfilters.hg19)
  library(svbams)
  library(Rsamtools)
  library(rtracklayer)
  library(graph)
  data(germline_filters)
  data(transcripts)
  ##
  ## read in some CNVs
  ##
  cv.extdata <- system.file("extdata", package="svcnvs")
  segs <- readRDS(file.path(cv.extdata, "cgov44t_segments.rds"))
  seqlevels(segs, pruning.mode="coarse") <- c("chr5", "chr8")
  extdata <- system.file("extdata", package="svbams")
  bview <- BamViews(bamPaths=file.path(extdata, "cgov44t_revised.bam"))
  ##
  ## Begin testing internals of sv_amplicons
  ##
  ## gr is the raw segmentation
  params <- ampliconParams()
  segs$is_amplicon <- segs$seg.mean > params[["AMP_THR"]]
  ag <- AmpliconGraph(ranges=segs,
                      filters=germline_filters,
                      params=params)
  if(FALSE){
    saveRDS(ag, file="AmpliconGraph0ada417.rds")
  }
  ag.100362c <- readRDS("AmpliconGraph100362c.rds")
  expect_identical(ag, ag.100362c)

  ag <- makeAGraph(segs, germline_filters, params)
  if(FALSE){
    saveRDS(ag, file="makeAGraphffab104.rds")
  }
})

