context("amplicons")

test_that("AmpliconGraph", {
  ag <- AmpliconGraph()
  expect_true(validObject(ag))
})

test_that("sv_amplicons_internals", {
  ## for system.file to work, devtools must be reloaded
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

  ##
  ## Begin testing internals of sv_amplicons
  ##
  amplicon_filters <- germline_filters
  params <- ampliconParams()
  ## gr is the raw segmentation
  segs <- gr

  af <- amplicon_filters
  af$border_size <- params$border_size
  af$overhang <- params$overhang
  AMP_THR <- params$AMP_THR
  ##AMP_THR <- amplicon_filters[["AMP_THR"]]
  segs$is_amplicon <- segs$seg.mean > AMP_THR
  ag <- AmpliconGraph(ranges=segs,
                      border_size=af[["border_size"]],
                      assembly_gaps=af[["assembly_gaps"]],
                      centromeres=af[["centromeres"]],
                      germline_cnv=af[["germline_cnv"]],
                      outliers=af[["outliers"]],
                      overhang=af[["overhang"]])
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

  centromeres <- af[["centromeres"]]
  ag <- trimRangesOverlappingCentromere(ag, centromeres)
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
  ##
  ## REFACTOR: could this step use the saved improper read pairs
  ##
  ## 
  rp <- svalignments::get_readpairs(ag, bamPaths(bview))
  ag1 <- addFocalDupsFlankingAmplicon(ag, rp, NULL)
  ag2 <- addFocalDupsFlankingAmplicon(ag, rp, LOW_THR)
  expect_identical(ag1, ag2)
  ag <- ag1

  qr1 <- focalAmpliconDupRanges(ag, LOW_THR=NULL, MAX_SIZE=500e3)
  qr2 <- focalAmpliconDupRanges(ag, LOW_THR=LOW_THR, MAX_SIZE=500e3)

  ##
  ## focalAmpliconDupRanges
  ##
  object <- ag
  g <- ampliconRanges(object)
  is_dup1 <- isDuplication(ranges(object), minimum_foldchange=NULL)
  is_dup2 <- isDuplication(ranges(object), minimum_foldchange=LOW_THR)
  ## is_dup1 is a length-0 logical vector
  expect_identical(is_dup1, logical())

  ## remove any lower copy ranges that are very large
  dup_g1 <- ranges(object)[is_dup1] ## is length-0 GRanges.  No ranges will be removed
  dup_g2 <- ranges(object)[is_dup2] ## identifies several large ranges
  dup_g <- dup_g[width(dup_g) < MAX_SIZE]
  g <- sort(c(dup_g, g))
  if(length(g) > 0){
    g <- expandGRanges(g, 1e3)
  }
  ##start(g) <- start(g)-1e3
  ##end(g) <- end(g)+1e3
  ## we do not want to link germline events
  germ <- germline(object)
  germ <- expandGRanges(germ, 10e3)
  ##start(germ) <- start(germ)-10e3
  ##end(germ) <- pmin(end(germ)+10e3, seqlengths(germ)[as.character(seqnames(germ))])
  germ <- reduce(germ)
  g <- filterBy(g, germ, type="within")
  reduce(g)



  queryRanges(ag1) <- qr1
  queryRanges(ag2) <- qr2

  irp <- get_improper_readpairs(ag, bamPaths(bview))
  ##
  ## At this point, focal duplications added to the graph have not
  ## been linked to any of the seeds
  ##
  ##paired_bin_filter <- af[["paired_bin_filter"]]
  ##param <- FilterEdgeParam(minimum_maxdist=50, bad_bins=paired_bin_filter)
  param <- FilterEdgeParam(minimum_maxdist=50, bad_bins=GRanges())
  ag <- linkFocalDups(ag, irp, LOW_THR=LOW_THR, edgeParam=param)
  ag <- linkAmplicons(ag, irp, edgeParam=param)

  ##
  ## burrough in on linkNearAmplicons
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
  ##
  ## end linkNearAmplicons


  ag <- linkNearAmplicons(ag, maxgap=500e3)
  expect_identical(ag, object)

  expected <- list("chr8:128,692,001",
                   c("chr5:176,034,001", "chr8:129,164,001", "chr8:129,353,001"),
                   c("chr8:129,353,001", "chr8:128,692,001"),
                   c("chr8:128,692,001", "chr8:129,164,001"))
  names(expected) <- c("chr5:176,034,001",
                       "chr8:128,692,001",
                       "chr8:129,164,001",
                       "chr8:129,353,001")
  expect_identical(e, expected)
  expect_identical(names(e), n)

  trace(sv_amplicons, browser)
  ag2 <- sv_amplicons(bview,
                     segs=gr,
                     amplicon_filters=germline_filters,
                     params=ampliconParams(),
                     transcripts=transcripts)
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
