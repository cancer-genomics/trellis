context("Revising sequence junctions")

.test_that <- function(name, expression) NULL

## overlapping hemizygous deletions on chromosome 3

.test_that("scratch", {
  library(svbams)
  bamdir <- system.file("extdata", package="svbams")
  bamfile <- file.path(bamdir, "cgov10t.bam")
  svfile <- file.path(bamdir, "cgov10t_deletions.rds")
  cbsfile <- file.path(bamdir, "cgov10t_cbs.rds")
  ##
  ## should regenerate from bamfile
  ##
  irpfile <- file.path(bamdir, "cgov10t_irp.rds")
  library(svfilters.hg19)
  data(germline_filters)
  filters <- reduceGenomeFilters(germline_filters)
  gr <- readRDS(cbsfile)
  gr <- gr[gr$seg.mean < log2(0.75)]

  file <- file.path("structuralvar/data/alignments/0improper", id.rds)
  aview <- AlignmentViews2(bview, file)
  file <- file.path("structuralvar/data/preprocess/3background_adj", id.rds)
  pview <- PreprocessViews2(bview)
  paths(pview) <- file
  setScale(pview) <- 1000
  set.seed(123)

  gr <- keepSeqlevels(gr, "chr3", pruning.mode="coarse")
  gr <- gr[3:5]
  names(gr) <- paste0("sv", 1:3)
  prp <- properReadPairs(bam_path=bamPaths(bview),
                         gr=gr, param=param)
  prp_index <- initializeProperIndex2(gr, prp, zoom.out=1)
  iparams <- improperAlignmentParams(which=gr, mapqFilter=0)
  irp <- getImproperAlignmentPairs(bamPaths(bview), iparams)
  names(irp) <- paste0("i", seq_along(irp))
  untrace(initializeImproperIndex2, browser)
  irp_index1 <- initializeImproperIndex2(cnv, irp, param)
  irp_index2 <- updateImproperIndex2(cnv, irp, maxgap=4e3)
  irp_index3 <- .match_index_variant(irp_index1, cnv, irp_index2)


  #irp <- addImproperReadPairs2(cnv, aview, param=param)


  irp_index1 <- initializeImproperIndex2(cnv, irp, param)
  irp_index2 <- updateImproperIndex2(cnv, irp, maxgap=2e3)
  irp_index3 <- .match_index_variant(irp_index1, cnv, irp_index2)

  sv <- sv_deletions(gr,
                     aview,
                     bview,
                     pview,
                     filters)
  saveRDS(sv, file=file.path("structuralvar/data/segment/1deletions", id.rds))

  sv.old <- readRDS(file.path("structuralvar/data/segment/1deletions", id.rds))
  ##sv.old <- sv.old[grep("chr3", as.character(seqnames(variant(sv.old))))]

  g.old <- variant(sv.old)
  g.old$seg.mean <- copynumber(sv.old)

  g.sv <- variant(sv)

  sv.old.diff <- sv.old[!overlapsAny(g.old, g.sv, type="equal")]
  g.delta <- variant(sv.old.diff)
  g.delta$seg.mean <- copynumber(sv.old.diff)
  g.delta$calls <- calls(sv.old.diff)

  ## chr3 region
  g.new <- subsetByOverlaps(variant(sv), g.delta[c("sv14", "sv15")])
  ## chrX region
  g.new <- subsetByOverlaps(variant(sv), g.delta[c("sv1", "sv58")])



  param <- DeletionParam()
  gr_filters <- filters
  cnv <- gr
  cnv <- germlineFilters(cnv, gr_filters, pview, param)
  thr <- log2(homozygousThr(param))
  cncalls <- ifelse(cnv$seg.mean < thr, "homozygous", "hemizygous")
  prp <- svalignments:::properReadPairs(bam_path=bamPaths(aview),
                                        gr=cnv, param=param)
  prp_index <- initializeProperIndex2(cnv, prp, zoom.out=1)
  irp <- addImproperReadPairs2(cnv, aview, param=param)
  irp_index1 <- initializeImproperIndex2(cnv, irp, param)
  irp_index2 <- updateImproperIndex2(cnv, irp, maxgap=2e3)
  irp_index3 <- .match_index_variant(irp_index1, cnv, irp_index2)

  sv <- deletion_call(aview, pview, gr, gr_filters)
  calls(sv) <- rpSupportedDeletions(sv, param=param, pview=pview)
  ##is_hemizygous <- calls(sv)=="hemizygous"
  ##sv <- sv[!is_hemizygous]
  sv <- reviseEachJunction(sv, pview, aview, param)
  sv2 <- removeSameStateOverlapping2(sv)
  sv <- sv2
  copynumber(sv) <- granges_copynumber(variant(sv), pview)
  calls(sv) <- rpSupportedDeletions(sv, param=param, pview=pview)
  if(FALSE){
    is_hemizygous <- calls(sv) == "hemizygous"
    sv <- sv[!is_hemizygous]
    if(length(sv) == 0) return(sv)xo
  }
  sv2 <- removeSameStateOverlapping2(sv)

  sv <- sv2
  indexImproper(sv) <- updateImproperIndex(sv, maxgap=500)
  calls(sv) <- rpSupportedDeletions(sv, param, pview=pview)
  is_hemizygous <- calls(sv) == "hemizygous"
  ##sv <- sv[!is_hemizygous]
  sv2 <- leftHemizygousHomolog(sv, pview, param)
  sv3 <- rightHemizygousHomolog(sv2, pview, param)
  calls(sv3) <- rpSupportedDeletions(sv3, param, pview)
  sv4 <- sv3[calls(sv3) != "hemizygous"]

  sv5 <- refineHomozygousBoundaryByHemizygousPlus(sv4)
  sv6 <- callOverlappingHemizygous(sv5)
  sv7 <- removeSameStateOverlapping2(sv6) 

  sv8 <- SVFilters(sv7, gr_filters, pview, param=param)
  sv9 <- groupSVs(sv8)
  id <- names(aview)
  sv9 <- allProperReadPairs(sv9, param, bfile=bamPaths(bview), zoom.out=1)

  chr3 <- keepSeqlevels(variant(sv), "chr3", pruning.mode="coarse")

  cncalls <- gsub("\\+", "", calls(sv))
  is.homdel <- cncalls == "homozygous"
  sv.homdel <- sv[is.homdel]
  sv.hemdel <- sv[!is.homdel]
  sv.hemdel <- adjudicateHemizygousOverlap2(sv.hemdel)
  sv.homdel <- adjudicateHomozygousOverlap2(sv.homdel)
})
