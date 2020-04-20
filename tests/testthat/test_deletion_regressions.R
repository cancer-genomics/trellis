context("Deletion regressions")

expect_identical2 <- function(sv1, sv2){
  variant(sv1) <- granges(variant(sv1))
  variant(sv2) <- granges(variant(sv2))
  expect_equivalent(sv1, sv2)
}

cgov44t_preprocess <- function(){
    data(bins1kb, package="svfilters.hg19", envir=environment())
  extdata <- system.file("extdata", package="svbams")
  id <- "cgov44t_revised.bam"
  bamfile <- file.path(extdata, id)
  path <- system.file("extdata", package="svbams")
  segs <- readRDS(file.path(path, "segs.4adcc78.rds"))

  ##gr <- readRDS("~/Dropbox/OvarianCellLines/structuralvar/data/segment/0cbs/CGOV44T.bam.rds")
  cnvpath <- system.file("extdata", package="svbams")
  gr <- readRDS(file.path(cnvpath, "cgov44t_segments.rds"))
  segs <- keepSeqlevels(gr, "chr15", pruning.mode="coarse")

  irp.file <- file.path(extdata, "cgov44t_improper.rds")
  irp <- readRDS(irp.file)
  ddir <- system.file("extdata", package="svbams",
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
  dp <- DeletionParam(remove_hemizygous=FALSE)
  dels <- sv_deletions(pdat)
  if(FALSE){
    saveRDS(dels, file="sv_deletions.ba3c739.rds")
  }
  path <- system.file("extdata", package="svbams")
  dels.ba3c739 <- readRDS(file.path(path, "sv_deletions.ba3c739.rds"))
  dels.ba3c739 <- rename(sort(dels.ba3c739))
  expect_equivalent(dels.ba3c739, dels)
})

test_that("deletion_call", {
  library(svfilters.hg19)
  pdat <- cgov44t_preprocess()
  result <- deletion_call(pdat)
  if(FALSE){
    saveRDS(result, file="deletion_call.4adcc78.rds")
  }
  path <- system.file("extdata", package="svbams")
  expected <- readRDS(file.path(path, "deletion_call.4adcc78.rds"))
  expect_equivalent(result, expected)
})

test_that("addImproperReadPairs2", {
  pdat <- cgov44t_preprocess()
  improper_rp <- pdat$read_pairs[["improper"]]
  mapq <- mcols(first(improper_rp))$mapq > 30 & mcols(last(improper_rp))$mapq > 30
  improper_rp <- improper_rp[mapq]
  ##
  ## TODO: this cutoff will be much too conservative in samples where tumor
  ## purity is less than 90%. Add tumor_purity to param object and take into
  ## account tumor_purity for determining cutoff
  ##
  cnv <- germlineFilters(pdat)
  irp <- improperRP(cnv, improper_rp)
  ##
  ## The current version return a GAlignmentPairs object with 2 fewer RPs. These
  ## additional RPs are more than 10kb from the candidate deletions -- they are
  ## excluded in "deletion_call.4adcc78.rds", so this irp object is correct even
  ## though it differs from irp.4adcc78
  if(FALSE){
    saveRDS(irp, file="addImproperReadPairs2.4adcc78.rds")
  }
  path <- system.file("extdata", package="svbams")
  irp.4adcc78 <- readRDS(file.path(path, "addImproperReadPairs2.4adcc78.rds"))
  expect_equivalent(irp, irp.4adcc78[1:72])
})

test_that("rpSupportedDeletions", {
  pdat <- cgov44t_preprocess()
  path <- system.file("extdata", package="svbams")
  sv <- readRDS(file.path(path, "deletion_call.4adcc78.rds"))
  calls <- rpSupportedDeletions(sv,
                                DeletionParam(),
                                pdat$bins)
  expect_identical(calls, "homozygous+")
})

test_that("rpSupportedDeletions_fails", {
  ## if the rpSupportedDeletion function filters variants that overlap with bins, then
  ## the calls returned
  pdat <- cgov44t_preprocess()
  path <- system.file("extdata", package="svbams")
  sv <- readRDS(file.path(path, "deletion_call.4adcc78.rds"))
  ##
  ## rpSupportedDeletions requires that each variant overlaps a bin. If one of
  ##  the variants does not overlap a bin, the calls vector returned will be
  ##  less than then legnth of the sv object
  ## -- below, reproduce the problem by removing all bins that overlap with the sv
  ##
  pdat$bins <- pdat$bins[!overlapsAny(pdat$bins, variant(sv))]
  calls <- rpSupportedDeletions(sv,
                                DeletionParam(),
                                pdat$bins)
  expect_error(stopifnot(length(calls) == length(sv)))
})




test_that("reviseEachJunction", {
  pdat <- cgov44t_preprocess()
  path <- system.file("extdata", package="svbams")
  sv <- readRDS(file.path(path, "deletion_call.4adcc78.rds"))
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
  path <- system.file("extdata", package="svbams")
  g.4adcc78 <- readRDS(file.path(path, "reviseEachJunction.4adcc78.rds"))
  expect_equivalent(g.4adcc78, g)
})

test_that("granges_copynumber", {
  pdat <- cgov44t_preprocess()
  path <- system.file("extdata", package="svbams")
  sv <- readRDS(file.path(path, "deletion_call.4adcc78.rds"))
  calls(sv) <- "homozygous+"
  path <- system.file("extdata", package="svbams")
  g <- readRDS(file.path(path, "reviseEachJunction.4adcc78.rds"))
  variant(sv) <- g
  cn <- granges_copynumber2(variant(sv), pdat$bins)
  expect_equal(-8.785, cn[[1]])
  copynumber(sv) <- cn
  expect_equal(copynumber(sv), cn)

  ## TODO maxgap should be part of the parameters
  elementMetadata(sv@improper)$names <- names(sv@improper)
  ##skip("Need help debugging GAlignments object")
  ## S4vectors no longer allows named GAlignments
  names(sv@improper)  <- NULL
  ##trace(updateImproperIndex, browser)
  index <- updateImproperIndex(sv, maxgap=500)
  if(FALSE){
    saveRDS(index, file="updateImproperIndex.4adcc78.rds")
    indexImproper(sv) <- index
    saveRDS(sv, file="sv_granges_copynumber.4adcc78.rds")
  }
  path <- system.file("extdata", package="svbams")
  index.4adcc78 <- readRDS(file.path(path, "updateImproperIndex.4adcc78.rds"))
  expect_equivalent(index, index.4adcc78)
})

## tests after granges_copynumber for a homozygous+ deletion
test_that("sv_deletions2", {
  library(svfilters.hg19)
  pdat <- cgov44t_preprocess()
  ## the only variant is homozygous+, so these functions are not doing anything
  path <- system.file("extdata", package="svbams")
  sv <- readRDS(file.path(path, "sv_granges_copynumber.4adcc78.rds"))
  sv3 <- finalize_deletions(sv, pdat)
  if(FALSE){
    saveRDS(sv3, file="allProperReadPairs.4adcc78.rds")
  }
  path <- system.file("extdata", package="svbams")
  sv.4adcc78 <- readRDS(file.path(path, "allProperReadPairs.4adcc78.rds"))
  expect_equivalent(sv.4adcc78, sv3)
})

test_that("germlineFilters", {
  library(svfilters.hg19)
  pdat <- cgov44t_preprocess()
  cnv <- germlineFilters(pdat)
  seqlevels(cnv, pruning.mode="coarse") <- "chr15"
  if(FALSE){
    saveRDS(cnvs, file="germlineFilters.9492f3f.rds")
  }
  path <- system.file("extdata", package="svbams")
  cnvs.9492f3f <- readRDS(file.path(path, "germlineFilters.9492f3f.rds"))
  expect_equivalent(cnv, cnvs.9492f3f)
})

.test_that <- function(name, expr) NULL

test_that("removeSameStateOverlapping2", {
  path <- system.file("extdata", package = "svbams")
  sv <- readRDS(file.path(path, "sv.rds"))
  skip("Need help debugging GAlignments object")
  sv2 <- removeSameStateOverlapping2(sv)
  expect_identical(length(sv2), 80L)
})

##
## Steps for confirming improper read pair alignment
##
## 1. extract the sequence of all improper read pairs in region near the borders
##    - for deletions with improperly paired reads (+ calls), borders can be narrowly defined
##
##    - for deletions without + calls, look in large region near segment boundary
##
## 2. write to fasta
## 3. blat alignment
## 4. check scores
## 5. update improper read pairs assigned to deletions object
## 6. update boundaries if needed
.test_that("blat_deletions", {
  full_bam <-  file.path("dcl01", "scharpf",
                         "data", "rscharpf",
                         "OvarianAlignments",
                         "hg19",
                         "eland",
                         "CGOV44T.bam")
  path <- system.file("extdata", package = "svbams")
  sv <- readRDS(file.path(path, "sv.rds"))
})


##
## Steps for evaluated split reads
##
## 1. look for mapped-unmapped near boundary
##
## 2. write sequences to fasta
##
## 3. blat alignment
##
## 4. check for split reads
##
## 5. revise boundaries if needed (++ call if split reads present?)
##
