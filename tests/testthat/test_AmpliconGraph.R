context("amplicons")

test_that("AmpliconGraph", {
  ag <- AmpliconGraph()
  expect_true(validObject(ag))
})

test_that("sv_amplicons_internals", {
  ## expectations defined from commit d2ed245
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
  gr <- readRDS(file.path(cv.extdata, "cgov44t_segments.rds"))
  extdata <- system.file("extdata", package="svbams")
  bview <- BamViews(bamPaths=file.path(extdata, "cgov44t_revised.bam"))
  amplicon_filters <- germline_filters
  params <- ampliconParams()
  ## gr is the raw segmentation
  segs <- gr
  ag2 <- sv_amplicons(bview,
                      segs=gr,
                      amplicon_filters=germline_filters,
                      params=ampliconParams(),
                      transcripts=transcripts)
  ##
  ## Begin testing internals of sv_amplicons
  ##
  ag <- makeAGraph(segs, amplicon_filters, params)
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
  tmp <- joinNearGRanges(ranges(ag), thr=0.05)
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
  ag <- addFocalDupsFlankingAmplicon(ag, rp, LOW_THR)
  qr <- focalAmpliconDupRanges(ag, LOW_THR=LOW_THR, MAX_SIZE=500e3)
  expected <- GRanges("chr5", IRanges(58519001, 58788001))
  seqinfo(expected) <- seqinfo(qr)
  expect_identical(qr[1], expected)
  queryRanges(ag) <- qr
  irp <- svalignments::get_improper_readpairs(ag, bamPaths(bview))
  ##
  ## At this point, focal duplications added to the graph have not
  ## been linked to any of the seeds
  ##
  param <- FilterEdgeParam(minimum_maxdist=50, bad_bins=GRanges())
  ag <- linkFocalDups(ag, irp, LOW_THR=LOW_THR, edgeParam=param)
  ag <- linkAmplicons(ag, irp, edgeParam=param)
  expect_identical(numEdges(ag), 2)
  expect_identical(numNodes(ag), 5L)

  ##
  ## linkNearAmplicons
  ##
  object <- ag
  maxgap <- 500e3
  ##
  ## *remove some arguments*
  ##
  hits <- findOverlaps(ampliconRanges(object), maxgap=maxgap)
  hits <- hits[!isSelfHit(hits)]
  hits <- hits[!isRedundantHit(hits)]
  if(length(hits)==0) return(object)
  new_edges <- paste(names(ampliconRanges(object))[queryHits(hits)],
                     names(ampliconRanges(object))[subjectHits(hits)],
                     sep="-")
  from <-  node1(new_edges)
  to <- node2(new_edges)
  existing <- edges(object)[from] ## is 'to' in any of the existing edges
  edge_exists <- mapply(function(to, existing) to %in% existing,
                        to=to, existing=existing)
  new_edges <- new_edges[!edge_exists]
  if(length(new_edges)==0) return(object)
  graph(object) <- addEdge(node1(new_edges),
                           node2(new_edges), graph(object))
  ## END

  ## TODO: add maxgap as parameter
  ag <- linkNearAmplicons(ag, maxgap=500e3)
  expect_identical(ag, object)
  e <-edges(ag)

  expected <- list("chr8:128,692,001",
                   c("chr8:129,353,001", "chr8:128,692,001", "chr8:129,164,001"),
                   c("chr5:176,034,001", "chr8:129,164,001", "chr8:129,353,001", "chr8:128,515,001"))
  names(expected) <- c("chr5:176,034,001",
                       "chr8:128,515,001",
                       "chr8:128,692,001")
  expect_identical(e[1:3], expected)
  expect_identical(nodes(ag), names(e))

  ag <- filterSmallAmplicons (ag)
  ag <- setAmpliconGroups (ag)
  ag <- setGenes (ag, transcripts)
  ag <- setDrivers (ag, transcripts, clin_sign=TRUE)
  ag <- setDrivers (ag, transcripts, clin_sign=FALSE)


  expect_identical(ag, ag2)
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
  library(rtracklayer)
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
