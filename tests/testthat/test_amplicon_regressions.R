context("Amplicon regressions")

test_that("AmpliconGraph", {
  ag <- AmpliconGraph()
  expect_true(validObject(ag))
})

test_that("sv_amplicons_internals", {
  ## expectations defined from commit d2ed245
  library(svfilters.hg19)
  library(svbams)
  library(Rsamtools)
  ##library(rtracklayer)
  library(graph)
  data(germline_filters)
  data(transcripts)
  ##
  ## read in some CNVs
  ##
  cv.extdata <- system.file("extdata", package="svcnvs")
  segs <- readRDS(file.path(cv.extdata, "cgov44t_segments.rds"))
  extdata <- system.file("extdata", package="svbams")
  bview <- BamViews(bamPaths=file.path(extdata, "cgov44t_revised.bam"))
  amplicon_filters <- germline_filters
  params <- ampliconParams()
  ##
  ## Begin testing internals of sv_amplicons
  ##
  ag <- makeAGraph(segs, amplicon_filters, params)
  ag.ffab104 <- readRDS("makeAGraphffab104.rds")
  expect_identical(ag, ag.ffab104)
  expect_true(validObject(ag))
  starts <- c(176034001, 128692001, 129164001, 129353001)
  ends <- c(177025001, 129163001, 129322001, 129615001)
  expected <- GRanges(c("chr5", rep("chr8", 3)),
                      IRanges(starts, ends))
  seqinfo(expected) <- seqinfo(transcripts)[c("chr5", "chr8"), ]
  isCircular(seqinfo(expected)) <- c(FALSE, FALSE)
  result <- granges(svclasses::amplicons(ag))
  names(result) <- NULL
  expect_identical(result, expected)
  expect_true(all(nodes(ag) %in% names(ranges(ag))))
  ##
  ## merge adjacent amplicons that have similar segment means
  ##
  tmp <- joinNearGRanges(ranges(ag), params)
  names(tmp) <- ampliconNames(tmp)

  ## the names of the nodes no longer correspond to the range names
  ## stopifnot(nodes(ag) %in% names(tmp)), and so
  ## setAmpliconGroups fails
  ranges(ag) <- tmp
  REMOTE <- file.exists(bamPaths(bview))
  if(!REMOTE) stop ("Path to BAM files is invalid")
  LOW_THR <- params[["LOW_THR"]]
  expect_true(LOW_THR > 0.75)
  ##
  ## REFACTOR: could this step use the saved improper read pairs
  ##
  ##
  rp <- svalignments::get_readpairs(ag, bamPaths(bview))
  ag <- addFocalDupsFlankingAmplicon(ag, rp, params)
  if(FALSE){
    saveRDS(ag, file="svcnvs/tests/testthat/addFocalDups.ffab104.rds")
  }
  ag.ffab104 <- readRDS("addFocalDups.ffab104.rds")
  expect_identical(ag, ag.ffab104)
})

test_that("focalAmpliconDupRanges", {
  params <- ampliconParams()
  ag <- readRDS("addFocalDups.ffab104.rds")
  ## TODO: ADD MAX_SIZE
  qr <- focalAmpliconDupRanges(ag, params)
  if(FALSE){
    saveRDS(qr, file="svcnvs/tests/testthat/focalAmpliconDup.ffab104.rds")
  }
  qr.ffab104 <- readRDS("focalAmpliconDup.ffab104.rds")
  expect_identical(qr, qr.ffab104)
  expected <- GRanges("chr5", IRanges(58519001, 58788001))
  seqinfo(expected) <- seqinfo(qr)
  expect_identical(qr[1], expected)
})

test_that("linkFocalDups", {
  ag <- readRDS("addFocalDups.ffab104.rds")
  qr <- readRDS("focalAmpliconDup.ffab104.rds")
  queryRanges(ag) <- qr
  extdata <- system.file("extdata", package="svbams")
  bview <- BamViews(bamPaths=file.path(extdata, "cgov44t_revised.bam"))
  irp <- svalignments::get_improper_readpairs(ag, bamPaths(bview))
  ##
  ## At this point, focal duplications added to the graph have not
  ## been linked to any of the seeds
  ##
  params <- ampliconParams()
  ag <- linkFocalDups(ag, irp, params)
  if(FALSE){
    saveRDS(ag, file="linkFocalDups.ffab104.rds")
  }
  ag.ffab104 <- readRDS("linkFocalDups.ffab104.rds")
  expect_identical(ag.ffab104, ag)
})

test_that("linkAmplicons", {
  ag <- readRDS("linkFocalDups.ffab104.rds")
  extdata <- system.file("extdata", package="svbams")
  bview <- BamViews(bamPaths=file.path(extdata, "cgov44t_revised.bam"))
  irp <- svalignments::get_improper_readpairs(ag, bamPaths(bview))
  params <- ampliconParams()
  ag <- linkAmplicons(ag, irp, edgeParam=params[["edge"]])
  if(FALSE){
    saveRDS(ag, file="linkAmplicons.4adcc78.rds")
  }
  ag.4adcc78 <- readRDS("linkAmplicons.4adcc78.rds")
  expect_identical(ag.4adcc78, ag)
})

test_that("linkNearAmplicons", {
  params <- ampliconParams()
  ag <- readRDS("linkAmplicons.4adcc78.rds")
  ag <- linkNearAmplicons(ag, maxgap=params[["maxgap"]])
  if(FALSE){
    saveRDS(ag, file="linkNearAmplicons.4adcc78.rds")
  }
  ag.4adcc78 <- readRDS("linkNearAmplicons.4adcc78.rds")
  expect_identical(ag.4adcc78, ag)
##  e <-edges(ag)
##  e <- e[-1]
##  expected <- list("chr8:128,692,001",
##                   c("chr8:129,353,001", "chr8:128,692,001", "chr8:129,164,001"),
##                   c("chr5:176,034,001", "chr8:129,164,001", "chr8:129,353,001",
##                     "chr8:128,515,001"))
##  names(expected) <- c("chr5:176,034,001",
##                       "chr8:128,515,001",
##                       "chr8:128,692,001")
##  expect_identical(e[1:3], expected)
##  expect_identical(nodes(ag), names(e))
})

test_that("filterSmallAmplicons", {
  ag <- readRDS("linkNearAmplicons.4adcc78.rds")
  ag <- filterSmallAmplicons (ag)
  if(FALSE){
    saveRDS(ag, file="filterSmallAmplicons.4adcc78.rds")
  }
  ag.4adcc78 <- readRDS("filterSmallAmplicons.4adcc78.rds")
  expect_identical(ag.4adcc78, ag)
})

test_that("setAmpliconGroups", {
  ag <- readRDS("filterSmallAmplicons.4adcc78.rds")
  ag <- setAmpliconGroups (ag)
  if(FALSE){
    saveRDS(ag, file="setAmpliconGroups.4adcc78.rds")
  }
  ag.4adcc78 <- readRDS("setAmpliconGroups.4adcc78.rds")
  expect_identical(ag.4adcc78, ag)
})

test_that("setGenes", {
  library(svfilters.hg19)
  ag <- readRDS("setAmpliconGroups.4adcc78.rds")
  ag <- setGenes (ag, transcripts)
  if(FALSE){
    saveRDS(ag, file="setAmpliconGenes.4adcc78.rds")
  }
  ag.4adcc78 <- readRDS("setAmpliconGenes.4adcc78.rds")
  expect_identical(ag.4adcc78, ag)
})

test_that("setDrivers", {
  library(svfilters.hg19)
  ag <- readRDS("setAmpliconGenes.4adcc78.rds")
  ag <- setDrivers (ag, transcripts, clin_sign=TRUE)
  ag <- setDrivers (ag, transcripts, clin_sign=FALSE)
  if(FALSE){
    saveRDS(ag, file="setDrivers.4adcc78.rds")
  }
  ag.4adcc78 <- readRDS("setDrivers.4adcc78.rds")
  expect_identical(ag.4adcc78, ag)
})

test_that("sv_amplicons", {
  library(svfilters.hg19)
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
  expect_identical(ag.4adcc78, ag2)
})

test_that("no germline filter", {
  library(svfilters.hg19)
  library(svbams)
  library(Rsamtools)
  library(graph)
  library(svalignments)
  data(germline_filters)
  data(transcripts)
  ##
  ## read in some CNVs
  ##
  cv.extdata <- system.file("extdata", package="svcnvs")
  segs <- readRDS(file.path(cv.extdata, "cgov44t_segments.rds"))
  extdata <- system.file("extdata", package="svbams")
  bview <- BamViews(bamPaths=file.path(extdata, "cgov44t_revised.bam"))
  params <- ampliconParams()
  germline_filters[["germline_cnv"]] <- GRanges()
  germline_filters[["outliers"]] <- GRanges()
  ##
  ## Begin testing internals of sv_amplicons
  ##
  ag <- makeAGraph(segs, germline_filters, params)
  merged <- joinNearGRanges(ranges(ag), params)
  names(merged) <- ampliconNames(merged)
  ranges(ag) <- merged
  rp <- get_readpairs(ag, bamPaths(bview))
  ag <- addFocalDupsFlankingAmplicon(ag, rp, params)
  queryRanges(ag) <- focalAmpliconDupRanges(ag, params)
  irp <- get_improper_readpairs(ag, bamPaths(bview))
  ag <- linkFocalDups(ag, irp, params)
  ag <- linkAmplicons(ag, irp, edgeParam=params[["edge"]])
  ag <- linkNearAmplicons(ag, maxgap=params[["maxgap"]])
  ag <- filterSmallAmplicons (ag)
  ag <- setAmpliconGroups (ag)
  ag <- setGenes (ag, transcripts)
  ag <- setDrivers (ag, transcripts, clin_sign=TRUE)
  ag <- setDrivers (ag, transcripts, clin_sign=FALSE)

  ag2 <- sv_amplicons(bview, segs,
                      germline_filters,
                      params, transcripts)
  expect_identical(ag, ag2)
  if(FALSE){
    saveRDS(ag2, file="sv_deletion.4adcc78.rds")
  }
  ag.4adcc78 <- readRDS("sv_deletion.4adcc78.rds")
  expect_identical(ag2, ag.4adcc78)
})

.test_that <- function(name, expr) NULL

.test_that("Full amplicon analysis of CGOV44T", {
  ##
  ## This test works for commmit 05d01c6
  ##
  ## saved amplicons graphs for the ovarian cell lines can be reproduced with
  ## commit 05d01c6
  library(svfilters.hg19)
  library(Rsamtools)
  ##library(rtracklayer)
  library(graph)
  library(svovarian)
  library(svcnvs)
  data(germline_filters)
  data(bviews_hg19)
  id <- "CGOV44T.bam"
  setwd("/dcl01/scharpf/data/rscharpf/projects/OvarianCellLines")
  gr <- readRDS("structuralvar/data/segment/0cbs/CGOV44T.bam.rds")
  bview <- bviews_hg19[, id]
  ag2 <- sv_amplicons(bview,
                      segs=gr,
                      amplicon_filters=germline_filters,
                      params=ampliconParams(),
                      transcripts=transcripts)

  ## this evaluates to TRUE prior for commmit 05d01c6 and earlier
  ## (LOW_THR is explicitly set to NULL)
  expected <- readRDS("structuralvar/data/segment/1amplicons/CGOV44T.bam.rds")
  expect_identical(ag2, expected)
})
