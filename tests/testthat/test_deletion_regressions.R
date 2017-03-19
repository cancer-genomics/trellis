context("Deletion regressions")

test_that("deletions_segs", {
  library(Rsamtools)
  library(svfilters.hg19)
  library(svpreprocess)
  data(germline_filters, package="svfilters.hg19")
  extdata <- system.file("extdata", package="svbams")
  id <- "CGOV44T.bam"
  id.rds <- paste0(id, ".rds")
  bamfile <- file.path(extdata, id)
  bview <- Rsamtools::BamViews(bamPaths=bamfile)
  data(segments, package="svcnvs")
  data(deletion, package="svcnvs")
  gr <- variant(deletion)
  gr <- expandGRanges(gr, 10000)
  segs <- keepSeqlevels(segments, "chr15", pruning.mode="coarse")
  segs <- segs[overlapsAny(segs, gr)]
  seqlevels(segs) <- seqlevels(gr)
  seqinfo(segs) <- seqinfo(gr)
  if(FALSE){
    saveRDS(segs, file="segs.4adcc78.rds")
  }
  segs.4adcc78 <- readRDS("segs.4adcc78.rds")
  expect_identical(segs, segs.4adcc78)

  ## seqinfo is required in the AlignmentViews object. Adding seqinfo to the
  ## BamViews ensures this information is propogated to the AlignmentViews
  ## object
  br <- bamRanges(bview)
  seqlevels(br) <- seqlevels(segs)
  seqinfo(br) <- seqinfo(segs)
  bamRanges(bview) <- br
  if(FALSE){
    saveRDS(bview, file="bview.4adcc78.rds")
  }
  if(FALSE){
    iparams <- improperAlignmentParams(which=gr, mapqFilter=30)
    irp <- getImproperAlignmentPairs(bview,
                                     iparams)
    ##
    ## TODO: refactor sv_deletions. Currently, we have to save the improperly
    ## paired reads to disk because sv_deletions reads from disk
    ##
    saveRDS(irp, file=irp.file)
  }
})

test_that("sv_deletions", {
  library(svfilters.hg19)
  bview <- readRDS("bview.4adcc78.rds")
  segs <- readRDS("segs.4adcc78.rds")
  ##
  ## extract improper alignments
  ##
  irp.file <- "getImproperAlignmentPairs.rds"
  aview <- AlignmentViews2(bview, path=irp.file)
  ddir <- system.file("extdata", package="svpreprocess",
                      mustWork=TRUE)
  cov.file <- file.path(ddir, "preprocessed_coverage.rds")
  data(pview, package="svpreprocess")
  paths(pview) <- cov.file

  filters <- reduceGenomeFilters(germline_filters)
  dparam <- DeletionParam()
  dels <- sv_deletions(gr=segs,
                       aview=aview,
                       bview=bview,
                       pview=pview,
                       gr_filters=filters,
                       param=dparam)
  if(FALSE){
    saveRDS(dels, file="sv_deletions.ba3c739.rds")
  }
  dels.ba3c739 <- readRDS("sv_deletions.ba3c739.rds")
  expect_identical(dels, dels.ba3c739)
})

test_that("deletion_call", {
  library(svfilters.hg19)
  filters <- reduceGenomeFilters(germline_filters)
  dparam <- DeletionParam()
  bview <- readRDS("bview.4adcc78.rds")
  segs <- readRDS("segs.4adcc78.rds")

  ddir <- system.file("extdata", package="svpreprocess",
                      mustWork=TRUE)
  cov.file <- file.path(ddir, "preprocessed_coverage.rds")
  data(pview, package="svpreprocess")
  paths(pview) <- cov.file

  ##
  ## deletion_call
  ##
  bview <- readRDS("bview.4adcc78.rds")
  segs <- readRDS("segs.4adcc78.rds")
  ##
  ## extract improper alignments
  ##
  irp.file <- "getImproperAlignmentPairs.rds"
  aview <- AlignmentViews2(bview, path=irp.file)
  result <- deletion_call(aview=aview,
                          pview=pview,
                          cnv=segs,
                          gr_filters=filters,
                          param=dparam)
  if(FALSE){
    saveRDS(result, file="deletion_call.4adcc78.rds")
  }
  expected <- readRDS("deletion_call.4adcc78.rds")
  expect_identical(result, expected)
})

test_that("addImproperReadPairs2", {
  library(svalignments)
  dparam <- DeletionParam()
  cnv <- readRDS("segs.4adcc78.rds")
  bview <- readRDS("bview.4adcc78.rds")
  irp.file <- "getImproperAlignmentPairs.rds"
  aview <- AlignmentViews2(bview, path=irp.file)
  ##
  ## TODO: this cutoff will be much too conservative in samples where tumor
  ## purity is less than 90%. Add tumor_purity to param object and take into
  ## account tumor_purity for determining cutoff
  ##
  irp <- addImproperReadPairs2(cnv, aview, param=dparam)
  if(FALSE){
    saveRDS(irp, file="addImproperReadPairs2.4adcc78.rds")
  }
  irp.4adcc78 <- readRDS("addImproperReadPairs2.4adcc78.rds")
  expect_identical(irp, irp.4adcc78)
})

getPview <- function(){
  ddir <- system.file("extdata", package="svpreprocess",
                      mustWork=TRUE)
  cov.file <- file.path(ddir, "preprocessed_coverage.rds")
  data(pview, package="svpreprocess")
  paths(pview) <- cov.file
  pview
}

test_that("rpSupportedDeletions", {
  sv <- readRDS("deletion_call.4adcc78.rds")
  pview <- getPview()
  dparam <- DeletionParam()
  calls <- rpSupportedDeletions(sv, param=dparam,
                                pview=pview)
  expect_identical(calls, "homozygous+")
})

test_that("reviseEachJunction", {
  sv <- readRDS("deletion_call.4adcc78.rds")  
  calls(sv) <- "homozygous+"
  pview <- getPview()
  bview <- readRDS("bview.4adcc78.rds")
  irp.file <- "getImproperAlignmentPairs.rds"
  aview <- AlignmentViews2(bview, path=irp.file)
  dparam <- DeletionParam()
  sv <- reviseEachJunction(object=sv,
                           pview=pview,
                           aview=aview,
                           param=dparam)
  sv <- removeSameStateOverlapping2(sv)
  g <- variant(sv)
  if(FALSE){
    saveRDS(g, file="reviseEachJunction.4adcc78.rds")
  }
  g.4adcc78 <- readRDS("reviseEachJunction.4adcc78.rds")
  expect_identical(g.4adcc78, g)
})

test_that("granges_copynumber", {
  pview <- getPview()
  sv <- readRDS("deletion_call.4adcc78.rds")  
  calls(sv) <- "homozygous+"
  g <- readRDS("reviseEachJunction.4adcc78.rds")
  variant(sv) <- g
  cn <- granges_copynumber(variant(sv), pview)
  expect_equal(-8.785, cn[[1]])
  copynumber(sv) <- cn
  expect_equal(copynumber(sv), cn)

  ## TODO maxgap should be part of the parameters
  index <- updateImproperIndex (sv, maxgap=500)
  if(FALSE){
    saveRDS(index, file="updateImproperIndex.4adcc78.rds")
    indexImproper(sv) <- index
    saveRDS(sv, file="sv_granges_copynumber.4adcc78.rds")
  }
  index.4adcc78 <- readRDS("updateImproperIndex.4adcc78.rds")
  expect_identical(index, index.4adcc78)
})

## tests after granges_copynumber for a homozygous+ deletion
test_that("sv_deletions2", {
  pview <- getPview()
  param <- DeletionParam()
  ## the only variant is homozygous+, so these functions are not doing anything
  sv <- readRDS("sv_granges_copynumber.4adcc78.rds")
  expect_identical(leftHemizygousHomolog(sv, pview, param), sv)
  expect_identical(rightHemizygousHomolog(sv, pview, param), sv)

  library(svfilters.hg19)
  filters <- reduceGenomeFilters(germline_filters)
  expect_identical(SVFilters(sv, filters, pview, param=param), sv)
  sv2 <- groupSVs(sv)
  expect_identical(groupedVariant(sv2), factor(1))

  bview <- readRDS("bview.4adcc78.rds")
  irp.file <- "getImproperAlignmentPairs.rds"
  aview <- AlignmentViews2(bview, path=irp.file)
  id <- names(aview)
  sv3 <- allProperReadPairs(sv2, param,
                            bfile=bamPaths(bview), zoom.out=1)
  prp <- proper(sv3)
  if(FALSE){
    saveRDS(sv3, file="allProperReadPairs.4adcc78.rds")
  }
  sv.4adcc78 <- readRDS("allProperReadPairs.4adcc78.rds")
  expect_identical(sv.4adcc78, sv3)
})

test_that("germlineFilters", {
  library(svfilters.hg19)
  dparam <- DeletionParam()
  filters <- reduceGenomeFilters(germline_filters)

  ddir <- system.file("extdata", package="svpreprocess",
                      mustWork=TRUE)
  cov.file <- file.path(ddir, "preprocessed_coverage.rds")
  data(pview, package="svpreprocess")
  paths(pview) <- cov.file

  data(segments, package="svcnvs")
  segs <- keepSeqlevels(segments, "chr15", pruning.mode="coarse")
  ##
  ## We only have a small region of the bam file
  ## - subset the segments
  data(deletion, package="svcnvs")
  gr <- variant(deletion)
  gr <- expandGRanges(gr, 10000)
  segs <- subsetByOverlaps(segs, gr)

  cnvs <- germlineFilters(cnv=segs,
                          germline_filters=filters,
                          pview=pview,
                          param=dparam)
  if(FALSE){
    saveRDS(cnvs, file="germlineFilters.9492f3f.rds")
  }
  cnvs.9492f3f <- readRDS("germlineFilters.9492f3f.rds")
  expect_identical(cnvs, cnvs.9492f3f)

  expected <- cnvs
  ##
  ## germlineFilters
  ##
  cnv <- segs
  not_germline <- isNotGermline(cnv, filters, dparam)
  ##
  ## idea is to compute the segment means for a much larger region and then to
  ##  calculate the fold change from the original cnv to the broader region -
  ##  for this unit test, the result should be 0
  egr <- expandGRanges(cnv, 0 * width(cnv))
  expect_identical(start(egr), start(cnv))
  cn <- granges_copynumber(egr, pview)
  expect_equal(cnv$seg.mean, as.numeric(cn),
               scale=1, tolerance=0.5)

  egr <- expandGRanges(cnv, 15 * width(cnv))
  cn_context <- granges_copynumber(egr, pview)
  fc <- (cnv$seg.mean - cn_context)
  expect_equivalent(fc < -3, c(FALSE, TRUE, FALSE))

  is_big <- isLargeHemizygous(cnv, dparam)
  expect_true(!any(is_big))
  select <- !is_big & not_germline & fc < -0.515
  cnv <- cnv[select]
  cnvr <- reduce(cnv)
  K <- width(cnvr) > 2000
  cnvr <- cnvr[K]
  cnv <- cnv[K]
  names(cnv) <- paste0("sv", seq_along(cnv))
  expect_identical(cnv, expected)
})

.test_that <- function(name, expr) NULL

test_that("removeSameStateOverlapping2", {
  sv <- readRDS("sv.rds")
  sv2 <- removeSameStateOverlapping2(sv)
  expect_identical(length(sv2), 80L)
})
