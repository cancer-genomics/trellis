context("get_read_pairs")

cgov44t_preprocess<- function(){
  extdata <- system.file("extdata", package="svbams")
  id <- "cgov44t_revised.bam"
  bamfile <- file.path(extdata, id)

  cv.extdata <- system.file("extdata", package="svcnvs")
  segments <- readRDS(file.path(cv.extdata, "cgov44t_segments.rds"))
  ##irp.file <- "getImproperAlignmentPairs.rds"
  ##aview <- AlignmentViews2(bview, path=irp.file)
  ##irp <- readRDS(irp.file)
  ddir <- system.file("extdata", package="svpreprocess",
                      mustWork=TRUE)
  lr <- readRDS(file.path(ddir, "preprocessed_coverage.rds"))/1000
  seqlevels(bins1kb, pruning.mode="coarse") <- paste0("chr", c(1:22, "X"))
  bins1kb$log_ratio <- lr

  irp.params <- improperAlignmentParams(mapqFilter=0)
  improper_rp <- getImproperAlignmentPairs(bamfile, irp.params)


  aparams <- ampliconParams()
  is.amp <- segments$seg.mean > aparams$AMP_THR
  amp.gr <- segments[is.amp]
  rp.params <- properAlignmentParams(which=reduce(amp.gr, min.gapwidth=2000),
                                     mapqFilter=0)
  proper_rp <- getProperAlignmentPairs(bamfile, rp.params)
  read_pairs <- list(proper=proper_rp,
                     improper=improper_rp)
  pdat <- preprocessData2(bam.file=bamfile,
                         genome=genome(segments)[[1]],
                         segments=segments,
                         read_pairs=read_pairs,
                         bins=bins1kb)
}

test_that("get_read_pairs",{
  ## input object must have queryRanges method defined
  library(Rsamtools)
  library(graph)
  library(svfilters.hg19)
  library(svbams)
  library(svalignments)
  pdata <- cgov44t_preprocess()
  test <- sv_amplicons2(pdata)

  cv.extdata <- system.file("extdata", package="svcnvs")
  segs <- readRDS(file.path(cv.extdata, "cgov44t_segments.rds"))
  ag <- readRDS("setDrivers.4adcc78.rds")
  extdata <- system.file("extdata", package="svbams")
  bview <- BamViews(bamPaths=file.path(extdata, "cgov44t_revised.bam"))
  params <- ampliconParams()
  ag2 <- sv_amplicons(bview=bview,
                      segs=segs,
                      amplicon_filters=germline_filters,
                      params=params,
                      transcripts=transcripts)
  ag.4adcc78 <- readRDS("setDrivers.4adcc78.rds")
  expect_identical(ag, ag2)
  expect_identical(test, ag)
})
