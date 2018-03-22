context("IKZF2 fusion")

.test_that <- function(nm, expr) NULL

test_that("sequenceJunctions", {
  library(trellis)
  library(BSgenome)
  extdata <- system.file("extdata", package="trellis")
  rlist <- readRDS(file.path(extdata, "rlist_5to3p.rds"))
  jxns <- readRDS(file.path(extdata, "rlist_jxns.rds"))
  coding_jxns <- codingJunctions(jxns, "hg19")
  expect_identical(genome(coding_jxns)[[1]], "hg19")
  expect_identical(length(coding_jxns), 44L)
  expect_true(all(coding_jxns$tx_name != ""))
  expect_true(all(coding_jxns$"3p"$tx_name != ""))
  extdata2 <- system.file("extdata", package="trellis")
  if(FALSE){
    saveRDS(coding_jxns, file=file.path(extdata2, "coding_jxns.rds"))
  }
  expected <- readRDS(file.path(extdata2, "coding_jxns.rds"))
  expect_equivalent(coding_jxns, expected)
  ikaros.rid <- "9136-9181"
  expect_true(ikaros.rid %in% coding_jxns$rid)
})

test_that("ikaros_fusion", {
  library(trellis)
  library(BSgenome)
  extdata <- system.file("extdata", package="trellis")
  rlist <- readRDS(file.path(extdata, "rlist_5to3p.rds"))
  extdata2 <- system.file("extdata", package="trellis")
  coding_jxns <- readRDS(file.path(extdata2, "coding_jxns.rds"))
  rid <- coding_jxns$rid
  rlist <- rlist[rid]

  tx.cds <- loadTxdbTranscripts(build)
  tx <- tx.cds[["transcripts"]]
  cds <- tx.cds[["cds"]]
  txdb <- tx.cds[["txdb"]]

  orgdb <- loadOrgDb()
  bs.pkg <- paste0("BSgenome.Hsapiens.UCSC.", build)
  genome <- getBSgenome(bs.pkg)
  ##loadGenomeData()

  ikaros.rid <- c("9136-9181", "9181-9136")
  rl <- rlist[ikaros.rid]
  ikaros.jxn <- coding_jxns[ikaros.rid]
  ##
  ## first orientation
  tx5p.id <- strsplit(ikaros.jxn$tx_name[1], ",")[[1]]
  tx3p.id <- strsplit(ikaros.jxn$"3p"$tx_name[1], ",")[[1]]
  cds.5p <- cds[tx5p.id]
  cds.3p <- cds[tx3p.id]

  ## sometimes different transcripts result in the exact same CDS after clipping
  tx5p.l2 <- lapply(cds.5p, clipFivePrime, jxn=ikaros.jxn[1])
  tx5p <- tx5p.l2[!duplicated(tx5p.l2)]
  tx3p.l2 <- lapply(cds.3p, clipThreePrime, jxn=ikaros.jxn$"3p"[1])
  expect_identical(duplicated(tx3p.l2), c(FALSE, FALSE))
  tx3p <- tx3p.l2[!duplicated(tx3p.l2)]
  grl <- expand.grid2(tx5p, tx3p)

  ##fusions <- fusion_CDS(r2, tx, cds)
  grl2 <- fuseCDS(rlist[ikaros.rid[1]],
                  jxn=coding_jxns[ikaros.rid[1]],
                  cds)
  expect_identical(grl, grl2[["fusion"]])
  if(FALSE){
    saveRDS(grl2, file=file.path(extdata2, "fuse_cds.rds"))
    saveRDS(tx5p, file=file.path(extdata2, "cds_5p_tumor.rds"))
  }
})
