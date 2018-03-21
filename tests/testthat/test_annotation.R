context("Annotating fusions")

data.raw <- function(){
  rlist <- readRDS("~/Dropbox/OvarianCellLines/structuralvar/data/rearrangements/3blat_unmapped/CGOV12T.bam.rds")
  r <- rlist[["1107-1108"]]
  saveRDS(r, file="cgov12t_rear.rds")
}


test_that("fusionTable", {
  r <- readRDS("cgov12t_rear.rds")
  build <- "hg19"
  library(org.Hs.eg.db)
  orgdb <- org.Hs.eg.db
  library(BSgenome.Hsapiens.UCSC.hg19)
  genome <- BSgenome.Hsapiens.UCSC.hg19
  library(TxDb.Hsapiens.UCSC.hg19.refGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.refGene
  tx <- transcripts(txdb)
  cds <- suppressWarnings(cdsBy(txdb, "tx", use.names=TRUE))
  tab <- fusionTable(robj=r,
                     txdb=txdb,
                     tx=tx,
                     cds=cds,
                     genome=genome,
                     orgdb=orgdb,
                     id="CGOV7T")
  expect_true(sum(is.na(tab$cdsStart.gene1)) == 6)
})


.test_that <- function(name, expr) NULL

.test_that("scratch", {
  path <- "~/Dropbox/OvarianCellLines/structuralvar/data/rearrangements"
  rlist <- readRDS(file.path(path, "3blat_unmapped/CGOV2T.bam.rds"))
  tab <- fusionList(rlist, id='cgov2t', build="hg19")
})
