.test_that <- function(nm, expr) NULL

test_that("fuseCDS_Rlist", {
  loadGenomeData()
  extdata <- system.file("extdata", package="svrearrange")
  rlist <- readRDS(file.path(extdata, "rlist_5to3p.rds"))
  extdata2 <- system.file("extdata", package="svfusions")
  coding_jxns <- readRDS(file.path(extdata2, "coding_jxns.rds"))
  rid <- coding_jxns$rid
  rlist <- rlist[rid]
  cds.list <- fuseCDS_Rlist(rlist, coding_jxns)
  fuse.list <- cds.list[["fusions"]]
  expect_identical(names(cds.list), c("fusions", "tum.5p",
                                      "tum.3p",
                                      "ref.5p", "ref.3p",
                                      "coding_junctions"))
  expect_identical(names(fuse.list), names(rlist))
  if(FALSE){
    saveRDS(cds.list, file=file.path(extdata2, "cds_list.rds"))
  }
  expect_is(fuse.list, "SimpleList")
  expect_is(fuse.list[[1]], "GRangesList")
  expect_identical(length(fuse.list), 44L)
})

## this takes too long
.test_that("translateCDS", {
  loadGenomeData()
  extdata <- system.file("extdata", package="svrearrange")
  extdata2 <- system.file("extdata", package="svfusions")
  rlist <- readRDS(file.path(extdata, "rlist_5to3p.rds"))
  coding_jxns <- readRDS(file.path(extdata2, "coding_jxns.rds"))
  coding_jxns <- coding_jxns[-grep("NR_", coding_jxns$tx_name)]
  rid <- coding_jxns$rid
  rlist <- rlist[rid]
  cds.list <- readRDS(file.path(extdata2, "cds_list.rds"))
  expect_true(all(elementNROWS(cds.list) == 43))
  if(FALSE){
    for(i in seq_along(cds.list)){
      cds.list[[i]] <- cds.list[[i]][1:10]
    }
  }
  fuse.list <- cds.list[["fusions"]]
  tum5p.list <- cds.list[["tum.5p"]]
  tum3p.list <- cds.list[["tum.3p"]]
  ref5p.list <- cds.list[["ref.5p"]]
  ref3p.list <- cds.list[["ref.3p"]]
  expect_identical(names(fuse.list), names(rlist))
  expect_identical(names(tum5p.list), names(rlist))
  fuse.p <- lapply(fuse.list, tumor_protein, genome=genome)
  if(FALSE){
    frame1 <- lapply(tprotein.list, "[[", 1)
    frame2 <- lapply(tprotein.list, "[[", 2)
    frame3 <- lapply(tprotein.list, "[[", 3)
  }
  tum5p.p <- lapply(tum5p.list, tumor_protein, genome=genome,
                    all.frames=TRUE)
  tum3p.p <- lapply(tum3p.list, tumor_protein, genome=genome,
                    all.frames=TRUE)
  ref5p.p <- lapply(ref5p.list, tumor_protein, genome=genome,
                    all.frames=TRUE)
  ref3p.p <- lapply(ref3p.list, tumor_protein, genome=genome,
                    all.frames=TRUE)
  proteins <- translateCDS(cds.list)
  expect_equivalent(proteins[["fusion"]], List(fuse.p))
  expect_equivalent(proteins[["tum5p"]], List(tum5p.p))  
  if(FALSE){
    saveRDS(proteins, file=file.path(extdata2, "cgov11t_proteins.rds"))
  }
  expect_identical(names(proteins[["fusion"]]), names(rlist))
  expect_identical(names(proteins[["tum5p"]]), names(coding_jxns)) 
})

test_that("pipeline", {
  loadGenomeData()
  library(svrearrange)
  extdata <- system.file("extdata", package="svrearrange")
  extdata2 <- system.file("extdata", package="svfusions")
  build <- "hg19"
  if(FALSE){
    rfile <- file.path(extdata, "CGOV11T_1.bam.rds")
    rlist <- readRDS(rfile)
    near.coding <- seqJunctionNearTx(rlist, build)
    rlist2 <- fiveTo3List(rlist[near.coding], build)
    ## not all jxns are coding
    jxns <- seqJunctions_Rlist(rlist2)
    coding_jxns <- codingJunctions(jxns, build)
    rlist <- rlist2[coding_jxns$rid]
    cds.list <- fuseCDS_Rlist(rlist, coding_jxns)

    tum.3p <- cds.list[["tum.3p"]][["8530-8533"]][[1]]
    ref.3p <- cds.list[["ref.3p"]][["8530-8533"]][[1]]
    jxn <- coding_jxns["8530-8533"]$"3p"
    expected <- clipThreePrime(ref.3p, jxn)
    expect_identical(tum.3p, expected)
    ##
    ## translate all CDS in three frames
    ##
    proteins <- translateCDS(cds.list)
  }
  rlist <- readRDS(file.path(extdata, "rlist_5to3p.rds"))
  extdata2 <- system.file("extdata", package="svfusions")
  coding_jxns <- readRDS(file.path(extdata2, "coding_jxns.rds"))
  rlist <- rlist[coding_jxns$rid]
  cds.list <- readRDS(file.path(extdata2, "cds_list.rds"))
  proteins <- readRDS(file.path(extdata2, "cgov11t_proteins.rds"))
  ##
  ## partition the amino acid sequence of the fusion into 5' and 3' components
  ## - repeat for each of the 3 possible reading frames
  ##
  partition.fusion <- partitionAASequence(proteins)
  expect_identical(names(partition.fusion), paste0("frame", 1:3))
  tum5p.frame1 <- lapply(proteins[["tum5p"]], "[[", 1)
  nAA.frame1 <- lapply(tum5p.frame1, width)
  frame1 <- lapply(proteins[["fusion"]], "[[", 1)
  expect_identical(elementNROWS(frame1),
                   elementNROWS(nAA.frame1)) 
  ##
  ## extract 3p portion from fused protein
  ##
  ref.frames <- organizeReferenceByFrame(proteins)
  expect_identical(names(ref.frames), paste0("frame", 1:3))
  ##
  ## evaluate whether in-frame by comparing to reference 3p protein
  ##
  inframe <- isInFrame_oneframe(fusion=partition.fusion[["frame1"]],
                                reference=ref.frames[["frame1"]])
  expect_equal(mean(unlist(inframe)), 0.33, tolerance=0.01)
  inframe.list <- inFrameList(fusion.frames=partition.fusion,
                              ref.frames=ref.frames)
  inframe <- inframe.list[[1]] | inframe.list[[2]] | inframe.list[[3]]

  ##trace(prematureStop_oneframe, browser)
  stops <- stopPositions_oneframe(partition.fusion[["frame1"]])
  nostops <- stops == "-1"
  inframe_and_nostops <- inframe & nostops
  ## there are 2 fusions that we would have considered inframe, but have premture stops
  ##table(unlist(inframe), unlist(nostops))
  fusions <- partition.fusion[["frame1"]]$fusion
  fusions.with.stops <- fusions[(inframe & !nostops)]
  fusions.with.stops <- fusions.with.stops[elementNROWS(fusions.with.stops) > 0]
  pos <- as.numeric(gregexpr("\\*", as.character(fusions.with.stops[[1]][[2]]))[[1]])
  ## the stop codon at 1032 is at the end of the 5-prime protein.
  ## - appropriate to filter this fusion
  expect_equivalent(pos, c(1032, 1982))
  expect_equal(mean(unlist(inframe_and_nostops)), 0.304, tolerance=0.01)
  nostop.list <- noPrematureStop(partition.fusion)
  number.valid <- sapply(lapply(nostop.list, unlist), sum)
  expect_equivalent(number.valid, c(24, 6, 7))
  ##
  ##  Its probably sufficient to always evaluate only the first frame
  ##
  expect_true(!any(unlist(!inframe.list[[1]] & inframe.list[[2]])))
  expect_true(!any(unlist(!inframe.list[[1]] & inframe.list[[3]])))
  ikaros.rid <- "9136-9181"
  expect_true(all(inframe.list[[1]][[ikaros.rid]]))
  inframe.nostop <- inFrameNoStop(nostop.list, inframe.list)
  ##trace(validFusions, browser)
  valid.fusions <- validFusions(partition.fusion,
                                cds.list,
                                inframe.nostop,
                                ref.frames)
  if(FALSE){
    extdata2 <- system.file("extdata", package="svfusions")
    saveRDS(valid.fusions, file=file.path(extdata2, "valid_fusions.rds"))
  }
  expected <- readRDS(file.path(extdata2, "valid_fusions.rds"))
  expect_equivalent(valid.fusions, expected)
  expect_true(ikaros.rid %in% names(valid.fusions[[1]]))
})

