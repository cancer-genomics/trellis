context("RearrangementList")

test_that("Rearrangement", {
  Rearrangement2 <- function(linkedBins=GRanges(),
                             improper=.GAlignmentPairs(),
                             partitioning=integer(),
                             link1id=character(),
                             link2id=character(),
                             tags=GRanges(),
                             tag_map_to_linked_bin=character(),
                             modal_rearrangement=character(),
                             percent_rearrangement=numeric(),
                             fraction_linking_tags=numeric(),
                             split_reads=GRanges()){
    new("Rearrangement",
        linkedBins=linkedBins,
        improper=improper,
        partitioning=partitioning,
        link1id=link1id,
        link2id=link2id,
        tags=tags,
        tag_map_to_linked_bin=tag_map_to_linked_bin,
        modal_rearrangement=modal_rearrangement,
        percent_rearrangement=percent_rearrangement,
        fraction_linking_tags=fraction_linking_tags,
        split_reads=GRanges())
  }
  r <- Rearrangement2()
  expect_true(validObject(r))
  expect_is(splitReads(r), "GRanges")
  splitReads(r) <- GRanges()
  expect_is(splitReads(r), "GRanges")
})

test_that("RearrangementList", {
  ## constructor
  rl <- RearrangementList()
  expect_true(validObject(rl))
  rl$something <- character()
  expect_identical(rl$something, character())
  expect_identical(nrow(colData(rl)), 0L)
  expect_identical(ncol(colData(rl)), 1L)
})

.test_that <- function(nm, expr) NULL
##why does this fail?
.test_that("rbind", {
  library(S4Vectors)
  df1 <- DataFrame(a=letters)
  df2 <- DataFrame(a=letters)
  result <- rbind(df1, df2)
  expect_is(result, "DataFrame")
})

test_that("c", {
  tmp <- c(RearrangementList(), RearrangementList())
  expect_is(tmp, "RearrangementList")
  data(rear_list)
  tmp <- c(rear_list[1:2], rear_list[3:4])
  expect_is(tmp, "RearrangementList")
  expect_equivalent(tmp, rear_list[1:4])
  ##
  ## NOTE: for some reason do.call(rbind, list_of_DataFrames) does not work
  ##  - had to coerce to data.frame and then back to DataFrame
  ##  - this causes the colData slot not to be identical
})

test_that("splitReads<-", {
  data(rear_list)
  nms <- names(rear_list)
  g <- replicate(4, GRanges())
  grl <- GRangesList(g)
  names(grl) <- head(nms, 4)
  splitReads(rear_list) <- grl

  ##
  ## accessor
  ##
  sr <- head(splitReads(rear_list), 4)
  expect_equivalent(sr, grl)
})


test_that("rlist from vignette",{
  extdata <- system.file("extdata", package="svbams")
  rfile <- file.path(extdata, "CGOV11T_1.bam.rds")
  rlist <- readRDS(rfile)
  build <- "hg19"
  near.coding <- seqJunctionNearTx(rlist, build)
  rlist <- rlist[ near.coding ]
  expect_true(validObject(rlist))
})


test_that("rlist for ovarian clearcell CGOV30T1", {
  extdata <- system.file("extdata", package="trellis")
  rlist <- readRDS(file.path(extdata, "rlist_cgov30t1.rds"))
  ##options(warn=2, error=recover)
  near.coding <- seqJunctionNearTx(rlist=rlist,
                                   build="hg19")
  expect_true(all(near.coding))
  rlist2 <- rlist[ near.coding ]
  if(length(rlist2) == 0){
    msg <- "No rearrangements near coding regions"
    if(interactive()) stop(msg) else stop(msg); q('no')
  }
  rlist2 <- fiveTo3List(rlist2, build="hg19")
  jxns <- seqJunctions_Rlist(rlist2)
  if(FALSE){
    labid <- strsplit(rds.file, "\\.") %>%
      sapply("[", 1)
    saveRDS(rlist2, file="~/Dropbox/Software/svpackages/trellis/inst/extdata/rlist_cgov30t1.rds")
  }
  coding_jxns <- codingJunctions(jxns, "hg19")
  rlist2 <- rlist2[ names(coding_jxns) ]
  cds.list <- fuseCDS_Rlist(rlist2, coding_jxns)
  proteins <- translateCDS(cds.list)
  partition.fusion <- partitionAASequence(proteins)
  nostop.list <- noPrematureStop(partition.fusion)
  ref.frames <- organizeReferenceByFrame(proteins)
  inframe.list <- inFrameList(fusion.frames=partition.fusion,
                              ref.frames=ref.frames)
  inframe.nostop <- inFrameNoStop(nostop.list, inframe.list)
  fusions <- validFusions(partition.fusion,
                          cds.list,
                          inframe.nostop,
                          ref.frames)
  tab <- fusionTable2(fusions) %>%
    as_tibble %>%
    filter(gene.5prime != gene.3prime)
})
