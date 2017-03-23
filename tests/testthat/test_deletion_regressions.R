context("Deletion regressions")

## test_that("deletions_segs", {
##   library(Rsamtools)
##   library(svfilters.hg19)
##   library(svpreprocess)
##   data(germline_filters, package="svfilters.hg19")
##   extdata <- system.file("extdata", package="svbams")
##   ##id <- "CGOV44T.bam"
##   id <- "cgov44t_revised.bam"
##   id.rds <- paste0(id, ".rds")
##   bamfile <- file.path(extdata, id)
##   bview <- Rsamtools::BamViews(bamPaths=bamfile)
##   data(segments, package="svcnvs")
##   data(deletion, package="svcnvs")
##   gr <- variant(deletion)
##   gr <- expandGRanges(gr, 10000)
##   segs <- keepSeqlevels(segments, "chr15", pruning.mode="coarse")
##   segs <- segs[overlapsAny(segs, gr)]
##   seqlevels(segs) <- seqlevels(gr)
##   seqinfo(segs) <- seqinfo(gr)
##   if(FALSE){
##     saveRDS(segs, file="segs.4adcc78.rds")
##   }
##   segs.4adcc78 <- readRDS("segs.4adcc78.rds")
##   expect_identical(segs, segs.4adcc78)
## 
##   ## seqinfo is required in the AlignmentViews object. Adding seqinfo to the
##   ## BamViews ensures this information is propogated to the AlignmentViews
##   ## object
##   br <- bamRanges(bview)
##   seqlevels(br) <- seqlevels(segs)
##   seqinfo(br) <- seqinfo(segs)
##   bamRanges(bview) <- br
##   if(FALSE){
##     saveRDS(bview, file="bview.4adcc78.rds")
##   }
## })



expect_identical2 <- function(sv1, sv2){
  variant(sv1) <- granges(variant(sv1))
  variant(sv2) <- granges(variant(sv2))
  expect_equivalent(sv1, sv2)
}

cgov44t_preprocess<- function(){
  extdata <- system.file("extdata", package="svbams")
  id <- "cgov44t_revised.bam"
  bamfile <- file.path(extdata, id)
  segs <- readRDS("segs.4adcc78.rds")
  irp.file <- "getImproperAlignmentPairs.rds"
  ##aview <- AlignmentViews2(bview, path=irp.file)
  irp <- readRDS(irp.file)
  ddir <- system.file("extdata", package="svpreprocess",
                      mustWork=TRUE)
  lr <- readRDS(file.path(ddir, "preprocessed_coverage.rds"))/1000
  seqlevels(bins1kb, pruning.mode="coarse") <- paste0("chr", c(1:22, "X"))
  bins1kb$log_ratio <- lr

  del.gr <- segs[segs$seg.mean < hemizygousThr(DeletionParam())]
  proper.del <- properReadPairs(bamfile,
                                gr=reduce(del.gr, min.gapwidth=2000))
  rps <- list(improper=irp, proper_del=proper.del)
  pdat <- preprocessData(bam.file=bamfile,
                         genome=genome(segs)[[1]],
                         segments=segs,
                         read_pairs=rps,
                         bins=bins1kb)
}

test_that("sv_deletions", {
  library(svfilters.hg19)
  pdat <- cgov44t_preprocess()
  dels <- sv_deletions(pdat)
  if(FALSE){
    saveRDS(dels, file="sv_deletions.ba3c739.rds")
  }
  dels.ba3c739 <- readRDS("sv_deletions.ba3c739.rds")
  dels.ba3c739 <- rename(sort(dels.ba3c739))
  expect_identical2(dels.ba3c739, dels)
})

test_that("deletion_call", {
  library(svfilters.hg19)
  pdat <- cgov44t_preprocess()
  result <- deletion_call(pdat)
  if(FALSE){
    saveRDS(result, file="deletion_call.4adcc78.rds")
  }
  expected <- readRDS("deletion_call.4adcc78.rds")
  expect_identical2(result, expected)
})

test_that("addImproperReadPairs2", {
  library(svalignments)
  pdat <- cgov44t_preprocess()
  rps <- pdat[["read_pairs"]]
  ##
  ## TODO: this cutoff will be much too conservative in samples where tumor
  ## purity is less than 90%. Add tumor_purity to param object and take into
  ## account tumor_purity for determining cutoff
  ##
  cnv <- germlineFilters(pdat)
  irp <- improperRP(cnv, rps$improper)
  ##
  ## The current version return a GAlignmentPairs object with 2 fewer RPs. These
  ## additional RPs are more than 10kb from the candidate deletions -- they are
  ## excluded in "deletion_call.4adcc78.rds", so this irp object is correct even
  ## though it differs from irp.4adcc78
  if(FALSE){
    saveRDS(irp, file="addImproperReadPairs2.4adcc78.rds")
  }
  irp.4adcc78 <- readRDS("addImproperReadPairs2.4adcc78.rds")
  expect_identical(irp, irp.4adcc78[1:72])
})

test_that("rpSupportedDeletions", {
  pdat <- cgov44t_preprocess()
  sv <- readRDS("deletion_call.4adcc78.rds")
  calls <- rpSupportedDeletions(sv,
                                DeletionParam(),
                                pdat$bins)
  expect_identical(calls, "homozygous+")
})

test_that("reviseEachJunction", {
  pdat <- cgov44t_preprocess()
  sv <- readRDS("deletion_call.4adcc78.rds")  
  calls(sv) <- "homozygous+"
  ## note, we can not use the improper read pairs stored in the sv object
  rps <- pdat$read_pairs
  irp <- rps$improper
  sv <- reviseEachJunction(sv,
                           pdat$bins,
                           irp)
  sv <- removeSameStateOverlapping2(sv)
  g <- variant(sv)
  if(FALSE){
    saveRDS(g, file="reviseEachJunction.4adcc78.rds")
  }
  g.4adcc78 <- readRDS("reviseEachJunction.4adcc78.rds")
  expect_identical(g.4adcc78, g)
})

test_that("granges_copynumber", {
  pdat <- cgov44t_preprocess()
  sv <- readRDS("deletion_call.4adcc78.rds")  
  calls(sv) <- "homozygous+"
  g <- readRDS("reviseEachJunction.4adcc78.rds")
  variant(sv) <- g
  cn <- granges_copynumber2(variant(sv), pdat$bins)
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
  library(svfilters.hg19)
  pdat <- cgov44t_preprocess()
  ## the only variant is homozygous+, so these functions are not doing anything
  sv <- readRDS("sv_granges_copynumber.4adcc78.rds")
  sv3 <- finalize_deletions(sv, pdat)
  if(FALSE){
    saveRDS(sv3, file="allProperReadPairs.4adcc78.rds")
  }
  sv.4adcc78 <- readRDS("allProperReadPairs.4adcc78.rds")
  expect_identical2(sv.4adcc78, sv3)
})

test_that("germlineFilters", {
  library(svfilters.hg19)
  pdat <- cgov44t_preprocess()
  cnv <- germlineFilters(pdat)
  seqlevels(cnv, pruning.mode="coarse") <- "chr15"
  if(FALSE){
    saveRDS(cnvs, file="germlineFilters.9492f3f.rds")
  }
  cnvs.9492f3f <- readRDS("germlineFilters.9492f3f.rds")
  expect_identical(cnv, cnvs.9492f3f)
})

.test_that <- function(name, expr) NULL

test_that("removeSameStateOverlapping2", {
  sv <- readRDS("sv.rds")
  sv2 <- removeSameStateOverlapping2(sv)
  expect_identical(length(sv2), 80L)
})
