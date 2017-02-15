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
  ##install_local("~/Software/svpackages/svbams")
  ##load_all("~/Software/svpackages/svbams")
  data(germline_filters)
  data(transcripts)
  ##
  ## read in some CNVs
  ##
  cv.extdata <- system.file("extdata", package="svcnvs")
  gr <- readRDS(file.path(cv.extdata, "cgov44t_segments.rds"))
  ##gr.bed <- import("../../../svbams/data-raw/cgov44t_revised.bed")
  ##gr.bed <- keepSeqlevels(gr.bed, c("chr5", "chr8"), pruning.mode="coarse")
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
  LOW_THR <- af[["LOW_THR"]]
  ##
  ## REFACTOR: could this step use the saved improper read pairs
  ##
  ## 
  rp <- svalignments::get_readpairs(ag, bamPaths(bview))
  ag <- addFocalDupsFlankingAmplicon(ag, rp, LOW_THR)
  queryRanges(ag) <- focalAmpliconDupRanges(ag, LOW_THR=LOW_THR, MAX_SIZE=500e3)
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


  ag <- filterSmallAmplicons (ag)
  ag <- setAmpliconGroups (ag)
  ag <- setGenes (ag, transcripts)
  ag <- setDrivers (ag, transcripts, clin_sign=TRUE)
  ag <- setDrivers (ag, transcripts, clin_sign=FALSE)
  e <- edges(ag)
  n <- nodes(ag)
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

  ag2 <- sv_amplicons(bview,
                     segs=gr,
                     amplicon_filters=germline_filters,
                     params=ampliconParams(),
                     transcripts=transcripts)
  expect_identical(ag, ag2)
})
