context("5 to 3-prime")

test_that("fiveTo3Prime", {
  extdata <- system.file("extdata", package="svbams")
  rfile <- file.path(extdata, "CGOV11T_1.bam.rds")
  rlist <- readRDS(rfile)
  r <- rlist[[1]]
  if(FALSE){
    for(i in seq_along(rlist)){
      r <- rlist[[i]]
      irp <- improper(r)
      first(irp) <- updateObject(irp@first)
      last(irp) <- updateObject(irp@last)
      improper(r) <- irp
      rlist[[i]] <- r
    }
    saveRDS(rlist, rfile)
  }
  r2 <- posNeg(r, "hg19")
  expect_identical(geneNames(r2), c("C3orf67-AS1,C3orf67", "FHIT"))
  r3 <- negPos(r, "hg19")
  df3 <- rearDataFrame(r3)
  if(FALSE) ggRearrange(df3)
  bins <- overlappingTranscripts(r, "hg19")
  orientations <- fiveTo3Prime(r, "hg19")
  id <- "10546-10582"
  r <- rlist[[id]]
  r2 <- fiveTo3Prime(r, "hg19")[[1]]
  df <- rearDataFrame(r2, "hg19")
  if(FALSE){
    ggRearrange(df)
  }
  ##id <- "12263-12576"
  r <- rlist[[1]]
  ##tx <- overlappingTranscripts(r, "hg19")
  ##expect_identical(as.character(tx$gene_name), c("CCDC58", "noncoding1"))
  ##df <- rearDataFrame(r, "hg19")
  ##expect_identical(levels(df$region), c("CCDC58", "noncoding1"))
  tx <- overlappingTranscripts(r, "hg19")
  is.valid <- check_splitreads(r)
  expect_identical(sum(!is.valid), 0L)
  splitReads(r) <- splitReads(r)[is.valid]
  ##df <- rearDataFrame(r, "hg19")
  if(FALSE)
    ggRearrange(df)
  tx <- overlappingTranscripts(rlist[[8]], "hg19")
  expect_identical(length(tx), 2L)

  ##IKZF2-ERBB4
  r <- rlist[["9136-9181"]]
  r2 <- fiveTo3Prime(r, "hg19")[[1]]
  r3 <- fiveTo3Prime(r, "hg19")[[2]]
  df <- rearDataFrame(r2, "hg19")
  df2 <- rearDataFrame(r3, "hg19")
  expect_identical(geneNames(r3),
                   c("IKZF2", "ERBB4,MIR548F2"))
  expect_identical(geneNames(r2),
                   rev(c("IKZF2", "ERBB4,MIR548F2")))
  if(FALSE){
    ggRearrange(df)
    ggRearrange(df2)
  }
})

test_that("noncoding", {
  extdata <- system.file("extdata", package="svbams")
  rfile <- file.path(extdata, "CGOV11T_1.bam.rds")
  rlist <- readRDS(rfile)
  r <- rlist[[5]]
  olist <- fiveTo3Prime(r, "hg19")
  df <- rearDataFrameList(olist)
  levs <- levels(df$region)
  expect_identical(levs[1], "5'-PDGFRB")
  expect_identical(levs[2], "3'-UBTF")
})

test_that("inversions", {
  extdata <- system.file("extdata", package="svbams")
  rfile <- file.path(extdata, "CGOV11T_1.bam.rds")
  rlist <- readRDS(rfile)
  ##
  ## inversions
  ##
  r <- rlist[["9258-33360"]]
  expect_true(isInversion(r))
  ri <- negativeInversion1(r, build="hg19")
  expect_identical(names(linkedBins(ri)), "9258-33360")
  rii <- negativeInversion2(r, build="hg19")
  expect_identical(names(linkedBins(rii)), "33360-9258")
  df <- rearDataFrame(ri, "hg19")
  df2 <- rearDataFrame(rii, "hg19")
  if(FALSE){
    ggRearrange(df2)
    ggRearrange(df)
  }
  o.list <- fiveTo3Prime(r, "hg19")
  expect_identical(names(o.list),
                   c("9258-33360", "33360-9258"))

  df <- rearDataFrameList(o.list)
  if(FALSE) ggRearrange(df)

  ##
  ## To do:  Every rearrangement should have 2 possible orientations
  ##   inversion + +
  ##
  r <- rlist[["9281-33369"]]
  ri <- positiveInversion1(r, build="hg19")
  expect_identical(geneNames(ri), c("USP37", "TUT1,MTA2"))
  df <- rearDataFrame(ri, "hg19")
  if(FALSE) ggRearrange(df)

  r2 <- swap_bin_order(ri)
  expect_identical(geneNames(r2), rev(c("USP37", "TUT1,MTA2")))

  rii <- positiveInversion2(r, build="hg19")
  expect_identical(geneNames(rii), rev(c("USP37", "TUT1,MTA2")))
  df <- rearDataFrame(rii, "hg19")
  if(FALSE) ggRearrange(df)
})

test_that("seqJunctionNearCoding", {
  library(trellis)
  library(BSgenome)
  extdata <- system.file("extdata", package="svbams")
  rfile <- file.path(extdata, "CGOV11T_1.bam.rds")
  rlist <- readRDS(rfile)
  ##rl <- rlist[ikaros.rid]
  near.coding <- seqJunctionNearTx(rlist, "hg19")
  expect_identical(sum(near.coding), 26L)
  if(FALSE){
    extdata <- system.file("extdata", package="svbams")
    saveRDS(near.coding, file=file.path(extdata, "near_coding.rds"))
  }
  ikaros.rid <- "9136-9181"
  expect_true(ikaros.rid %in% names(rlist[near.coding]))
})

test_that("five_to_three", {
  library(BSgenome)
  extdata <- system.file("extdata", package="svbams")
  rfile <- file.path(extdata, "CGOV11T_1.bam.rds")
  rlist <- readRDS(rfile)
  rlist2 <- fiveTo3List(rlist, build="hg19")
  expect_equal(length(names(rlist2)), length(rlist2))
  if(FALSE){
    saveRDS(rlist2, file=file.path(extdata, "rlist_5to3p.rds"))
  }
  expect_is(rlist2, "RearrangementList")
  expect_equal(length(rlist2), 52)
  expect_is(rlist2[1], "RearrangementList")
  ikaros.rid <- "9136-9181"
  expect_true(ikaros.rid %in% names(rlist2))
})

test_that("seqJunctions_Rlist", {
  library(trellis)
  library(BSgenome)
  extdata <- system.file("extdata", package="svbams")
  rlist <- readRDS(file.path(extdata, "rlist_5to3p.rds"))
  jxns <- seqJunctions_Rlist(rlist)
  expect_is(jxns, "GRanges")
  if(FALSE){
    saveRDS(jxns, file=file.path(extdata, "rlist_jxns.rds"))
  }
})

.test_that <- function(name, expr)  NULL
.test_that("fiveTo3List", {
  rfile <- "~/Dropbox/OvarianCellLines/rearrangements/3blat_unmapped/CGOV1T.bam.rds"
  build <- "hg19"
  rlist <- readRDS(rfile)
  is.complex <- isComplex(rlist)
  rlist <- rlist[ !is.complex ]
  n.sr <- numberSplitReads(rlist)
  rlist <- rlist[ n.sr >= 1 ]
  near.coding <- seqJunctionNearTx(rlist, build)
  rlist <- rlist[ near.coding ]

  r <- rlist[["14704-14706"]]
  r2 <- RearrangementList(fiveTo3Prime(r, build))

  trace(rearDataFrame, browser)
  df <- rearDataFrame(r2[[1]], build, 5000)
  df <- rearDataFrame(r2[[2]], build, 5000)
  trace(rearDataFrameList, browser)
  df <- rearDataFrameList(r2)

  ggRearrange(df)

  linked <- linkedTo(r2)
  trace(rearDataFrameList, browser)
  df <- rearDataFrameList(r2)
  if(FALSE){
    ggRearrange(df)
  }
  rlist2 <- fiveTo3List(rlist, build="hg19")
})
