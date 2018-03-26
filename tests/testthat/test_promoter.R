context("Fusion involving promoter of 5' gene")

promoterFusion <- function(){
  ##
  ## Put the rearrangement junction in the promoter of YAP1
  ##
  library(GenomicRanges)
  library(svfilters.hg19)
  library(trellis)
  library(TxDb.Hsapiens.UCSC.hg19.refGene)
  data(transcripts, package="svfilters.hg19", envir=environment())
  tx <- get("transcripts")
  is.yap1 <- mcols(tx)$gene_name == "YAP1"
  yap1 <- reduce(tx[is.yap1])
  tx.yap1 <- subsetByOverlaps(tx, yap1)[1]
  yap1.promoter <- promoters(tx.yap1, upstream=5000, downstream=0)
  txdb <- TxDb.Hsapiens.UCSC.hg19.refGene
  tx <- transcripts(txdb)

  data(rear_cgov7t, package="trellis")
  ##rlist <- rear_cgov7t
  ##attributes(class(rlist))$package <- "trellis"
  ##for(i in seq_along(rlist)){
  ##  r <- rlist@data[[i]]
  ##  attributes(class(r))$package <- "trellis"
  ##  rlist@data[[i]] <- r
  ##}
  ##rear_cgov7t <- rlist
  ##save(rear_cgov7t, file="data/rear_cgov7t.rda")
  rear.id <- "16978-17458"
  rear <- rear_cgov7t[[rear.id]]

  irp <- improper(rear)
  i1 <- subjectHits(findOverlaps(yap1.promoter, first(irp), maxgap=130000))

  F <- first(irp)[i1]
  delta <- start(F) - (start(yap1.promoter) + 3000)
  F2 <- F
  F2@start <- as.integer(start(F) - delta)
  all(overlapsAny(F2, yap1.promoter))
  first(irp)[i1] <- F2

  library(GenomicAlignments)
  i2 <- subjectHits(findOverlaps(yap1.promoter, last(irp), maxgap=130000))
  L <- last(irp)[i2]
  delta <- start(L) - (start(yap1.promoter) + 3000)
  L2 <- L
  L2@start <- as.integer(start(L) - delta)
  all(overlapsAny(L2, yap1.promoter))
  irp@last[i2] <- L2
  rear@improper <- irp
  ##
  lb <- linkedBins(rear)
  lb$linked.to <- reduce(as(L2, "GRanges"))
  linkedBins(rear) <- lb
  rear
}

test_that("promoter", {
  library(trellis)
  rear <- promoterFusion()
  txdb <- loadTxDb("hg19")
  tx <- transcripts(txdb)
  bs.pkg <- paste0("BSgenome.Hsapiens.UCSC.", "hg19")
  genome <- getBSgenome(bs.pkg)
  orgdb <- loadOrgDb()
  cds <- suppressWarnings(cdsBy(txdb, "tx", use.names=TRUE))

  trx <- rearrangementTranscripts(linkedBins(rear), tx, cds)
  identifier(trx) <- names(linkedBins(rear))
  left.cluster <- granges(linkedBins(rear))
  right.cluster <- granges(linkedBins(rear)$linked.to)
  modal.call <- modalRearrangement(rear)
  cdfun <- getCDFun(modal.call)
  cds.obj <- cdfun(left.cluster, right.cluster, trx)
  cds.fusions <- getCDS(rear, tx, cds)
  expect_identical(cds.obj, cds.fusions)
  expect_true(isFusion(cds.fusions))
  expect_identical(length(txC(cds.fusions)), 9L)

  clipped <- clip(cds.fusions)
  fused <- fuse(clipped)
  expect_identical(names(fused), "NM_001130145(promoter)::NM_032427")
  ## there are no CDS in YAP1
  fused.proteins <- tumorProtein(genome, fused)
  expect_true(all(fused$orientation=="right"))

  tx.nms <- unique(unlist(strsplit(names(fused.proteins), "::")))
  cds <- fullTranscripts(cds.fusions)
  ref.proteins <- referenceProtein(genome, cds, tx.nms)
  in_frame <- inFrameFusions(fused.proteins, ref.proteins, fused)
  expect_true(in_frame)
  tab <- .fusionTable(fused.txlist=fused,
                      fused.proteins=fused.proteins,
                      in_frame=in_frame,
                      org.db=orgdb,
                      txdb=txdb,
                      linkedbin.id=names(rear))
  expect_identical(tab$gene1, "YAP1")
  expect_identical(tab$chr.gene1, "chr11")
  fused.tx <- fusionList(RearrangementList(rear))
  expect_identical(tab, fused.tx)
  expect_true(fused.tx$inframe)
})

test_that("promoter2", {
  ## The second gene is not necessarily in-frame
  data(rear_cgov7t, package="trellis")
  rid <- "16350-16406"
  tab <- fusionList(rear_cgov7t[rid])
  expect_true(!tab$inframe)
})
