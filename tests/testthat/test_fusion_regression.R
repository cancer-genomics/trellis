context("Regression tests")

.test_that <- function(nm, expr) NULL

test_that("overlapsAnyTranscript", {
  data(rear_cgov7t, package="trellis")
  overlaps.tx <- overlapsAnyTranscript(rear_cgov7t, "hg19", maxgap=5000)
  expect_true(sum(overlaps.tx) == 24)
  if(FALSE){
    saveRDS(overlaps.tx, file="overlaps.tx.rds")
  }
  extdata <- system.file("extdata", package="svbams")
  expected <- readRDS(file.path(extdata, "overlaps.tx.rds"))
  expect_identical(overlaps.tx, expected)
})

.test_that("fusionList", {
  data(rear_cgov7t, package="trellis")
  test <- fusionList(rear_cgov7t, "CGOV7T")
  if(FALSE){
    saveRDS(test, file="fusionList.641eb81.rds")
  }
  extdata <- system.file("extdata", package="svbams")
  expected <- readRDS(file.path(extdata, "fusionList.641eb81.rds"))
  expect_identical(test, expected)
})

test_that("fusionTable", {
  library(org.Hs.eg.db)
  orgdb <- org.Hs.eg.db
  library(BSgenome.Hsapiens.UCSC.hg19)
  genome <- BSgenome.Hsapiens.UCSC.hg19
  library(TxDb.Hsapiens.UCSC.hg19.refGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.refGene
  data(rear_cgov7t)
  tx <- transcripts(txdb)
  cds <- suppressWarnings(cdsBy(txdb, "tx", use.names=TRUE))
  rid <- "16978-17458"
  tab <- fusionTable(robj=rear_cgov7t[[rid]],
                     txdb=txdb,
                     tx=tx,
                     cds=cds,
                     genome=genome,
                     orgdb=orgdb,
                     id="CGOV7T")
  extdata <- system.file("extdata", package="svbams")
  expected <- readRDS(file.path(extdata, "fusionList.641eb81.rds"))
  expected <- expected[expected$rearrangement.id == rid, ]
  expect_equivalent(tab, expected)
})

.test_that("getCds", {
  library(TxDb.Hsapiens.UCSC.hg19.refGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.refGene
  data(rear_cgov7t)
  rid <- "16978-17458"
  r <- rear_cgov7t[[rid]]
  tx <- transcripts(txdb)
  cds <- suppressWarnings(cdsBy(txdb, "tx", use.names=TRUE))
  cds.fusions <- getCDS(r, tx, cds)
  if(FALSE){
    saveRDS(cds.fusions, file="cds.fusions.677b50c.rds")
  }
  extdata <- system.file("extdata", package="svbams")
  expected <- readRDS(file.path(extdata, "cds.fusions.677b50c.rds"))
  expect_equivalent(cds.fusions, expected)
})

.test_that("fuse_and_clip", {
  ##
  ## a promoter is fused to two possible transcripts on chr11
  ##
  cds.fusions <- readRDS("cds.fusions.677b50c.rds")
  clipped <- clip(cds.fusions)
  fused <- fuse(clipped)
  if(FALSE){
    saveRDS(fused, file="fused.677b50c.rds")
  }
  expected <- readRDS("fused.677b50c.rds")
  expect_identical(fused, expected)
})


.test_that("tumorProtein", {
  library(BSgenome.Hsapiens.UCSC.hg19)
  bs_genome <- BSgenome.Hsapiens.UCSC.hg19
  fused.txlist <- readRDS("fused.677b50c.rds")
  fused.proteins <- tumorProtein(bs_genome, fused.txlist)
  if(FALSE){
    saveRDS(fused.proteins, file="fused.proteins.677b50c.rds")
  }
  expected <- readRDS("fused.proteins.677b50c.rds")
  expect_equivalent(fused.proteins, expected)
})

.test_that("fullTranscripts", {
  cds.fusions <- readRDS("cds.fusions.677b50c.rds")
  cds <- fullTranscripts(cds.fusions)
  if(FALSE){
    saveRDS(cds, file="cds.full.tx.677b50c.rds")
  }
  expected <- readRDS("cds.full.tx.677b50c.rds")
  expect_identical(cds, expected)
})

.test_that("referenceProtein", {
  library(BSgenome.Hsapiens.UCSC.hg19)
  bs_genome <- BSgenome.Hsapiens.UCSC.hg19
  cds <- readRDS("cds.full.tx.677b50c.rds")
  fused.proteins <- readRDS("fused.proteins.677b50c.rds")
  tx.nms <- unique(unlist(strsplit(names(fused.proteins), "::")))
  ref.protein <- referenceProtein(bs_genome, cds, tx.nms)
  if(FALSE){
    saveRDS(ref.protein, file="ref.protein.677b50c.rds")
  }
  expected <- readRDS("ref.protein.677b50c.rds")
  expect_equivalent(ref.protein, expected)
})

.test_that("inFrameFusions", {
  library(org.Hs.eg.db)
  library(BSgenome.Hsapiens.UCSC.hg19)
  ref.protein <- readRDS("ref.protein.677b50c.rds")
  fused.txlist <- readRDS("fused.677b50c.rds")
  fused.proteins <- readRDS("fused.proteins.677b50c.rds")
  in_frame <- inFrameFusions(fused.proteins, ref.protein,
                             fused.txlist)
  expect_true(all(in_frame))
})

.test_that(".fusionTable", {
  library(org.Hs.eg.db)
  orgdb <- org.Hs.eg.db
  library(BSgenome.Hsapiens.UCSC.hg19)
  bs_genome <- BSgenome.Hsapiens.UCSC.hg19
  library(TxDb.Hsapiens.UCSC.hg19.refGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.refGene
  fused.txlist <- readRDS("fused.677b50c.rds")
  fused.proteins <- readRDS("fused.proteins.677b50c.rds")
  data(rear_cgov7t)
  rid <- "16978-17458"
  robj <- rear_cgov7t[[rid]]
  tab <- .fusionTable(fused.txlist=fused.txlist,
                      fused.proteins=fused.proteins,
                      in_frame=rep(TRUE, 9),
                      org.db=orgdb,
                      txdb=txdb,
                      linkedbin.id=names(robj),
                      id="CGOV7T")
  expected <- readRDS("fusionList.641eb81.rds")
  expected <- expected[expected$rearrangement.id == rid, ]
  expect_equivalent(tab, expected)
})
