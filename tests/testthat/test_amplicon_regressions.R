context("Amplicon regressions")

test_that("AmpliconGraph", {
  ag <- AmpliconGraph()
  expect_true(validObject(ag))
})



cgov44t_preprocess <- function(){
  extdata <- system.file("extdata", package="svbams")
  id <- "cgov44t_revised.bam"
  bamfile <- file.path(extdata, id)

  irp.file <- file.path(extdata, "cgov44t_improper.rds")
  irp <- readRDS(irp.file)
  ddir <- system.file("extdata", package="trellis",
                      mustWork=TRUE)
  lr <- readRDS(file.path(ddir, "preprocessed_coverage.rds"))/1000
  seqlevels(bins1kb, pruning.mode="coarse") <- paste0("chr", c(1:22, "X"))
  bins1kb$log_ratio <- lr

  path <- system.file("extdata", package="trellis")
  segs <- readRDS(file.path(path, "cgov44t_segments.rds"))
  seqlevels(segs, pruning.mode="coarse") <- c("chr5", "chr8")
  ##amp.gr <- segs[segs$seg.mean < ampliconParams()$AMP_THR]
  ##proper.amp <- properReadPairs(bamfile,
  ##                              gr=reduce(amp.gr, min.gapwidth=2000))
  rps <- list(improper=irp)
  pdat <- preprocessData(bam.file=bamfile,
                         genome=genome(segs)[[1]],
                         segments=segs,
                         read_pairs=rps,
                         bins=bins1kb)
}

test_that("sv_amplicons", {
  library(Rsamtools)
  library(svfilters.hg19)
  ##
  ##  standard setup
  ##
  cv.extdata <- system.file("extdata", package="trellis")
  segs <- readRDS(file.path(cv.extdata, "cgov44t_segments.rds"))
  seqlevels(segs, pruning.mode="coarse") <- c("chr5", "chr8")
  extdata <- system.file("extdata", package="svbams")
  bview <- BamViews(bamPaths=file.path(extdata, "cgov44t_revised.bam"))
  params <- ampliconParams()
  tx <- loadTx("hg19")
  ##
  ## call sv_amplicons with a bunch of arguments
  ##
  ag2 <- sv_amplicons(bview=bview,
                      segs=segs,
                      amplicon_filters=germline_filters,
                      params=params,
                      transcripts=tx)

  path <- system.file("extdata", package="trellis")
  ag.4adcc78 <- readRDS(file.path(path, "setDrivers.4adcc78.rds"))

  expect_identical(ag.4adcc78, ag2)
  ##
  ## proposed setup.  Call sv_amplicons2 with single argument
  ##
  pdat <- cgov44t_preprocess()
  expect_identical(pdat$segments, segs)
  ag3 <- sv_amplicons2(pdat)
  ##
  ## these will not be exactly identical because the set of improper read pairs
  ## changes slightly.  However, the amplicon graph is exactly the same
  ##
  expect_identical(ampliconRanges(ag2), ampliconRanges(ag3))
  expect_identical(queryRanges(ag2), queryRanges(ag3))
  expect_identical(edges(ag2), edges(ag3))
})

test_that("initialize_graph", {
  ##
  ## standard
  ##
  cv.extdata <- system.file("extdata", package="trellis")
  segs <- readRDS(file.path(cv.extdata, "cgov44t_segments.rds"))
  seqlevels(segs, pruning.mode="coarse") <- c("chr5", "chr8")
  path <- system.file("extdata", package = "trellis")
  ag <- readRDS(file.path(path, "setDrivers.4adcc78.rds"))
  extdata <- system.file("extdata", package="svbams")
  bview <- BamViews(bamPaths=file.path(extdata, "cgov44t_revised.bam"))
  params <- ampliconParams()
  amplicon_filters <- germline_filters
  params <- ampliconParams()
  ag <- initialize_graph(segs, amplicon_filters, params)
  if(FALSE){
    saveRDS(ag, file="initialize_graph.a4d7744.rds")
  }
  path <- system.file("extdata", package = "trellis")
  ag.a4d7744 <- readRDS(file.path(path, "initialize_graph.a4d7744.rds"))
  expect_identical(ag, ag.a4d7744)
  ## proposed
  pdat <- cgov44t_preprocess()
  pdat$segments <- amplified_segments(pdat$segments, params)
  ag2 <- initialize_graph2(pdat, ampliconFilters(pdat$genome),
                           ampliconParams())
  expect_identical(ag, ag2)
})

test_that("add_amplicons", {
  path <- system.file("extdata", package = "trellis")
  ag <- readRDS(file.path(path, "addFocalDups.ffab104.rds"))
  path <- system.file("extdata", package = "trellis")
  query.ranges <- readRDS(file.path(path, "focalAmpliconDupRanges.a4d7744.rds"))
  queryRanges(ag) <- query.ranges
  expected <- ag
  ## proposed
  pdat <- cgov44t_preprocess()
  path <- system.file("extdata", package = "trellis")
  ag <- readRDS(file.path(path, "initialize_graph.a4d7744.rds"))
  ag2 <- add_amplicons(ag, pdat$bam.file, ampliconParams())
  if(FALSE){
    saveRDS(ag2, file="add_amplicons.a4d7744.rds")
  }
  expect_identical(ag2, expected)
})


test_that("link_amplicons", {
  ##
  ## standard
  ##
  extdata <- system.file("extdata", package="svbams")
  bam.file <- file.path(extdata, "cgov44t_revised.bam")
  path <- system.file("extdata", package = "trellis")
  ag <- readRDS(file.path(path, "add_amplicons.a4d7744.rds"))
  irp <- get_improper_readpairs(ag, bam.file)
  params <- ampliconParams()
  ag1 <- link_amplicons(ag, irp, params)
  if(FALSE){
    saveRDS(ag1, file="link_amplicons.a4d744.rds")
  }
  path <- system.file("extdata", package = "trellis")
  ag.a4d744 <- readRDS(file.path(path, "link_amplicons.a4d744.rds"))
  expect_identical(ag1, ag.a4d744)
  ##
  ## proposed
  ##
  pdat <- cgov44t_preprocess()
  improper_rp2 <- pdat$read_pairs[["improper"]]
  ag2 <- link_amplicons(ag, improper_rp2, params)
  expect_identical(ampliconRanges(ag2), ampliconRanges(ag1))
  expect_identical(queryRanges(ag2), queryRanges(ag1))
  expect_identical(edges(ag2), edges(ag1))
  ##
  ## if we use all the improper read pairs with a mapq filter we recover one fewer edge
  ##
  iparams <- improperAlignmentParams(mapqFilter=30,
    what=c("flag", "mrnm", "mpos"))
  improper_rp3 <- getImproperAlignmentPairs(pdat$bam.file,
                                            param=iparams, 
                                            build="hg19")
  ag3 <- link_amplicons(ag, improper_rp3, params)
  expect_lt(numEdges(ag3), numEdges(ag2))
})

test_that("annotate_amplicons", {
  path <- system.file("extdata", package = "trellis")
  expected <- readRDS(file.path(path, "setDrivers.4adcc78.rds"))
  path <- system.file("extdata", package = "trellis")
  ag <- readRDS(file.path(path, "link_amplicons.a4d744.rds"))
  tx <- loadTx("hg19")
  ag <- annotate_amplicons(ag, tx)
  expect_identical(ag, expected)
})


## ag3 is not identical to ag2
## -- check the internal functions
test_that("makeAGraph", {
  library(Rsamtools)
  library(svfilters.hg19)
  ##
  ##  standard setup
  ##
  cv.extdata <- system.file("extdata", package="trellis")
  segs <- readRDS(file.path(cv.extdata, "cgov44t_segments.rds"))
  seqlevels(segs, pruning.mode="coarse") <- c("chr5", "chr8")
  path <- system.file("extdata", package = "trellis")
  ag <- readRDS(file.path(path, "setDrivers.4adcc78.rds"))
  extdata <- system.file("extdata", package="svbams")
  bview <- BamViews(bamPaths=file.path(extdata, "cgov44t_revised.bam"))
  params <- ampliconParams()
  amplicon_filters <- germline_filters
  params <- ampliconParams()
  ##
  ## Begin testing internals of sv_amplicons
  ##
  ag <- makeAGraph(segs, amplicon_filters, params)
  path <- system.file("extdata", package = "trellis")
  ag.ffab104 <- readRDS(file.path(path, "makeAGraphffab104.rds"))
  expect_identical(ag, ag.ffab104)
  expect_true(validObject(ag))
  ##
  ## Proposed setup
  ##
  pdat <- cgov44t_preprocess()
  ## generates an error
  segs <- pdat$segments
  segs$is_amplicon <- segs$seg.mean > params$AMP_THR
  agnew <- makeAGraph2(segs, amplicon_filters, params)
  expect_identical(agnew, ag)
})

## the proposed and standard setup are the same here
test_that("joinNearGRanges", {
  ##
  ##  standard setup
  ##
  path <- system.file("extdata", package = "trellis")
  ag <- readRDS(file.path(path, "makeAGraphffab104.rds"))
  merged <- joinNearGRanges(ranges(ag), ampliconParams())
  if(FALSE){
    saveRDS(merged, file="merged.a4d7744.rds")
  }
  path <- system.file("extdata", package = "trellis")
  merged.a4d7744 <- readRDS(file.path(path, "merged.a4d7744.rds"))
  expect_identical(merged, merged.a4d7744)
})

test_that("get_readpairs", {
  ##
  ## standard
  ##
  path <- system.file("extdata", package = "trellis")
  ag <- readRDS(file.path(path, "initialize_graph.a4d7744.rds"))
  extdata <- system.file("extdata", package="svbams")
  bam.file <- file.path(extdata, "cgov44t_revised.bam")
  rp <- get_readpairs(ag, bam.file)
  ##
  ## proposed
  ##
  rp2 <- get_readpairs2(queryRanges(ag), bam.file)
  expect_identical(rp, rp2)
})

test_that("addFocalDups", {
  path <- system.file("extdata", package = "trellis")
  ag <- readRDS(file.path(path, "initialize_graph.a4d7744.rds"))
  ## standard and proposed are the same
  extdata <- system.file("extdata", package="svbams")
  bam.file <- file.path(extdata, "cgov44t_revised.bam")
  rp <- get_readpairs(ag, bam.file)
  ag <- addFocalDupsFlankingAmplicon(ag, rp, ampliconParams())
  path <- system.file("extdata", package = "trellis")
  ag.ffab104 <- readRDS(file.path(path, "addFocalDups.ffab104.rds"))
  expect_identical(ag, ag.ffab104)
  query.ranges <- focalAmpliconDupRanges(ag, ampliconParams())
  if(FALSE){
    saveRDS(query.ranges, file="focalAmpliconDupRanges.a4d7744.rds")
  }
  path <- system.file("extdata", package = "trellis")
  query.a4d7744 <- readRDS(file.path(path, "focalAmpliconDupRanges.a4d7744.rds"))
  expect_identical(query.ranges, query.a4d7744)
})


test_that("linkNearAmplicons", {
  params <- ampliconParams()
  path <- system.file("extdata", package = "trellis")
  ag <- readRDS(file.path(path, "linkAmplicons.4adcc78.rds"))
  ag <- linkNearAmplicons(ag, maxgap=params[["maxgap"]])
  if(FALSE){
    saveRDS(ag, file="linkNearAmplicons.4adcc78.rds")
  }
  path <- system.file("extdata", package = "trellis")
  ag.4adcc78 <- readRDS(file.path(path, "linkNearAmplicons.4adcc78.rds"))
  expect_identical(ag.4adcc78, ag)
})

test_that("filterSmallAmplicons", {
  path <- system.file("extdata", package = "trellis")
  ag <- readRDS(file.path(path, "linkNearAmplicons.4adcc78.rds"))
  ag <- filterSmallAmplicons (ag)
  if(FALSE){
    saveRDS(ag, file="filterSmallAmplicons.4adcc78.rds")
  }
  path <- system.file("extdata", package = "trellis")
  ag.4adcc78 <- readRDS(file.path(path, "filterSmallAmplicons.4adcc78.rds"))
  expect_identical(ag.4adcc78, ag)
})

test_that("setAmpliconGroups", {
  path <- system.file("extdata", package = "trellis")
  ag <- readRDS(file.path(path, "filterSmallAmplicons.4adcc78.rds"))
  ag <- setAmpliconGroups (ag)
  if(FALSE){
    saveRDS(ag, file="setAmpliconGroups.4adcc78.rds")
  }
  path <- system.file("extdata", package = "trellis")
  ag.4adcc78 <- readRDS(file.path(path, "setAmpliconGroups.4adcc78.rds"))
  expect_identical(ag.4adcc78, ag)
})

test_that("setGenes", {
  library(svfilters.hg19)
  path <- system.file("extdata", package = "trellis")
  ag <- readRDS(file.path(path, "setAmpliconGroups.4adcc78.rds"))
  tx <- loadTx("hg19")
  ag <- setGenes (ag, tx)
  if(FALSE){
    saveRDS(ag, file="setAmpliconGenes.4adcc78.rds")
  }
  path <- system.file("extdata", package = "trellis")
  ag.4adcc78 <- readRDS(file.path(path, "setAmpliconGenes.4adcc78.rds"))
  expect_identical(ag.4adcc78, ag)
})

test_that("setDrivers", {
  library(svfilters.hg19)
  path <- system.file("extdata", package = "trellis")
  ag <- readRDS(file.path(path, "setAmpliconGenes.4adcc78.rds"))
  tx <- loadTx("hg19")
  ag <- setDrivers (ag, tx, clin_sign=TRUE)
  ag <- setDrivers (ag, tx, clin_sign=FALSE)
  if(FALSE){
    saveRDS(ag, file="setDrivers.4adcc78.rds")
  }
  path <- system.file("extdata", package = "trellis")
  ag.4adcc78 <- readRDS(file.path(path, "setDrivers.4adcc78.rds"))
  expect_identical(ag.4adcc78, ag)
})

test_that("no germline filter", {
  library(svfilters.hg19)
  library(svbams)
  library(Rsamtools)
  library(graph)
  data(germline_filters)
  data(transcripts)
  ##
  ## read in some CNVs
  ##
  cv.extdata <- system.file("extdata", package="trellis")
  segs <- readRDS(file.path(cv.extdata, "cgov44t_segments.rds"))
  seqlevels(segs, pruning.mode="coarse") <- c("chr5", "chr8")
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
  tx <- loadTx("hg19")
  ag <- setGenes (ag, tx)
  ag <- setDrivers (ag, tx, clin_sign=TRUE)
  ag <- setDrivers (ag, tx, clin_sign=FALSE)

  ag2 <- sv_amplicons(bview, segs,
                      germline_filters,
                      params, tx)
  expect_identical(ag, ag2)
  if(FALSE){
    saveRDS(ag2, file="sv_deletion.4adcc78.rds")
  }
  path <- system.file("extdata", package = "trellis")
  ag.4adcc78 <- readRDS(file.path(path, "sv_deletion.4adcc78.rds"))
  expect_identical(ag2, ag.4adcc78)
})
