context("Deletion regressions")

test_that("sv_deletions", {
  library(Rsamtools)
  library(svfilters.hg19)
  library(svpreprocess)
  data(germline_filters, package="svfilters.hg19")
  extdata <- system.file("extdata", package="svbams")
  id <- "CGOV44T.bam"
  id.rds <- paste0(id, ".rds")
  bamfile <- file.path(extdata, id)
  bview <- Rsamtools::BamViews(bamPaths=bamfile)
  br <- bamRanges(bview)

  data(segments, package="svcnvs")
  segs <- keepSeqlevels(segments, "chr15", pruning.mode="coarse")

  ## seqinfo is required in the AlignmentViews object. Adding seqinfo to the
  ## BamViews ensures this information is propogated to the AlignmentViews
  ## object
  seqlevels(br) <- seqlevels(segs)
  seqinfo(br) <- seqinfo(segs)
  bamRanges(bview) <- br



  data(deletion, package="svcnvs")
  gr <- variant(deletion)
  gr <- expandGRanges(gr, 10000)
  segs <- segs[overlapsAny(segs, gr)]
  ##
  ## extract improper alignments
  ##
  iparams <- improperAlignmentParams(which=gr)
  irp <- getImproperAlignmentPairs(bview,
                                   iparams,
                                   mapq_thr=30,
                                   use.mcols=TRUE)
  ##
  ## TODO: refactor sv_deletions. Currently, we have to save the improperly
  ## paired reads to disk because sv_deletions reads from disk
  ##
  irp.file <- "getImproperAlignmentPairs.rds"
  saveRDS(irp, file=irp.file)
  aview <- AlignmentViews2(bview, path=irp.file)
  ##
  ## Create pviews object
  ##
  ddir <- system.file("extdata", package="svpreprocess",
                      mustWork=TRUE)
  cov.file <- file.path(ddir, "preprocessed_coverage.rds")
  data(pview, package="svpreprocess")
  paths(pview) <- cov.file

  filters <- reduceGenomeFilters(germline_filters)
  dparam <- DeletionParam()
  seqlevels(segs) <- seqlevels(gr)
  seqinfo(segs) <- seqinfo(gr)
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

  ##
  ## sv_deletions
  ##
  gr_filters <- filters
  gr <- segs
  param <- dparam

  ##
  ## deletion_call
  ##
  cnv <- germlineFilters(gr, filters, pview, dparam)
  ##
  ## TODO: this cutoff will be much too conservative in samples where tumor
  ## purity is less than 90%. Add tumor_purity to param object and take into
  ## account tumor_purity for determining cutoff
  ##
  thr <- log2(homozygousThr(dparam))
  cncalls <- ifelse(cnv$seg.mean < thr, "homozygous", "hemizygous")
  prp <- properReadPairs(bam_path=bamPaths(aview),
                         gr=cnv, param=param)
  prp_index <- initializeProperIndex2(cnv, prp, zoom.out=1)
  irp <- addImproperReadPairs2(cnv, aview, param=dparam)
  irp_index1 <- initializeImproperIndex2(cnv, irp, dparam)
  irp_index2 <- updateImproperIndex2(cnv, irp, maxgap=2e3)
  irp_index3 <- .match_index_variant(irp_index1, cnv, irp_index2)
  ##irp <- irp[validPairForDeletion(irp)]
  ##irp_index <- initializeImproperIndex2(cnv, irp, param)
  sv <- StructuralVariant(variant=cnv,
                          proper=prp,
                          improper=irp,
                          copynumber=cnv$seg.mean,
                          calls=cncalls,
                          index_proper=prp_index,
                          index_improper=irp_index3)

  sv2 <- deletion_call(aview, pview, gr, gr_filters)
  expect_identical(sv, sv2)

  calls(sv) <- rpSupportedDeletions(sv, param=dparam, pview=pview)
  is_hemizygous <- calls(sv)=="hemizygous"
  sv <- sv[!is_hemizygous]
  if(length(sv) == 0) return(sv)
  sv <- reviseEachJunction(sv, pview, aview, dparam)
  if(length(sv) == 0) return(sv)
  copynumber(sv) <- granges_copynumber(variant(sv), pview)
  calls(sv) <- rpSupportedDeletions(sv, param=param, pview=pview)
  expect_identical(calls(sv), "homozygous+")

  is_hemizygous <- calls (sv) == "hemizygous"
  sv <- sv[!is_hemizygous]
  if(length(sv) == 0) return(sv)
  sv <- removeSameStateOverlapping(sv)

  indexImproper(sv) <- updateImproperIndex (sv, maxgap=500)
  calls(sv) <- rpSupportedDeletions(sv, param, pview=pview)
  is_hemizygous <- calls(sv) == "hemizygous"
  sv <- sv[!is_hemizygous]
  if(length(sv) == 0) return(sv)

  sv2 <- leftHemizygousHomolog(sv, pview, param)
  sv3 <- rightHemizygousHomolog(sv2, pview, param)
  calls(sv3) <- rpSupportedDeletions(sv3, param, pview)
  message("Removing hemizygous deletions without rearranged RPs")
  sv4 <- sv3[calls(sv3) != "hemizygous"]
  if(length(sv4)==0) return(sv4)

  message("Refining homozygous boundaries by spanning hemizygous+")
  sv5 <- refineHomozygousBoundaryByHemizygousPlus(sv4)
  sv6 <- callOverlappingHemizygous(sv5)
  sv7 <- removeSameStateOverlapping(sv6) 

  sv8 <- SVFilters(sv7, gr_filters, pview, param=param)
  sv9 <- groupSVs(sv8)
  id <- names(aview)
  sv9 <- allProperReadPairs(sv9, param, bfile=bamPaths(bview), zoom.out=1)
  if(length(sv9@proper) > 25e3){
    proper(sv9) <- sv9@proper[sample(seq_along(sv9@proper), 25e3)]
    indexProper(sv9) <- initializeProperIndex3(sv9, zoom.out=1)
  }
  expect_identical(sv9, dels.ba3c739)
})
