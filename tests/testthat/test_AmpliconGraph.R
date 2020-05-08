context("AmpliconGraph")

test_that("AmpliconGraph", {
  ag <- AmpliconGraph()
  expect_true(validObject(ag))

  library(svfilters.hg19)
  library(svbams)
  library(Rsamtools)
  library(rtracklayer)
  library(graph)
  data(germline_filters, package="svfilters.hg19")
  data(transcripts, package="svfilters.hg19")
  ##
  ## read in some CNVs
  ##
  cv.extdata <- system.file("extdata", package="svbams")
  segs <- readRDS(file.path(cv.extdata, "cgov44t_segments.rds"))
  seqlevels(segs, pruning.mode="coarse") <- c("chr5", "chr8")
  extdata <- system.file("extdata", package="svbams")
  bview <- BamViews(bamPaths=file.path(extdata, "cgov44t_revised.bam"))
  ##bview <- BamViews(bamPaths=file.path(extdata, "cgov44t_test.bam"))
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
  path <- system.file("extdata", package="svbams")
  ag.100362c <- readRDS(file.path(path, "AmpliconGraph100362c.rds"))
  expect_equivalent(ag, ag.100362c)

  ag <- makeAGraph(segs, germline_filters, params)
  if(FALSE){
    saveRDS(ag, file="makeAGraphffab104.rds")
  }
})

