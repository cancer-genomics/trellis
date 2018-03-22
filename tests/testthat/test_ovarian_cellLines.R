context("Ovarian cell line fusions")

.test_that <- function(name, expr) NULL

##
## fusion in CGOV44T
test_that("cgov44t", {
  library(trellis)
  library(svalignments)
  data(rlist, package="trellis")
  test <- fusionList(rlist, id="CGOV44T")
  ## check that we've recovered the correct orientation'
  gene1.coords <- GRanges(test$chr.gene1, IRanges(test$cdsStart.gene1,
                                                  test$cdsEnd.gene1))
  ## the first gene should overlap the "linked.to" field
  expect_true(all(overlapsAny(gene1.coords, linkedTo(rlist)[2], maxgap=10000)))
  expect_true(!any(overlapsAny(gene1.coords, linkedBins(rlist)[1], maxgap=10000)))

  if(FALSE){
    extdata <- system.file("extdata", package="svalignments")
    unmap.file <- file.path(extdata, "blat_unmapped.txt")
    blat_unmap <- readBlat(unmap.file)
    split_reads <- rearrangedReads(rlist, blat_unmap, 500)
    ## put rearrangement in correct orientation
    df.list <- organizeReads(rlist, split_reads)
    dat <- df.list[[2]]
    dat$type2 <- paste0(dat$strand, dat$read)
    dat$side <- ifelse(dat$side=="left", "right somatic", "left somatic")
    dat$side <- factor(dat$side, levels=c("left somatic", "right somatic"))
    left.jxn <- max(dat$end[dat$type=="split read" & dat$side=="left somatic"])
    right.jxn <- min(dat$start[dat$type=="split read" & dat$side=="right somatic"])
    df <- data.frame(junction=c(left.jxn, right.jxn),
                     side=c("left somatic", "right somatic"))
    library(ggplot2)
    ggplot(dat, aes(xmin=start, xmax=end,
                      ymax=rpid+0.25, ymin=rpid-0.25)) +
      geom_rect(alpha=0.5, aes(color=type2, fill=type2)) +
      geom_vline(data=df, aes(xintercept=junction), linetype="dashed",
                 inherit.aes=FALSE) +
      theme(axis.text.x=element_text(size=7)) +
      facet_grid(~side, scales="free_x")
  }
})

##
## fusion in CGOV7T
##
test_that("Check YAP1-MAMML",{
  ## YAP1-MAMML
  data(rear_cgov7t)
  ##trace(fusionList, browser)
  fusions <- fusionList(rlist=rear_cgov7t,
                        id="CGOV7T")
  yap1_mamml <- fusions[fusions$gene1=="YAP1", ]
  expect_true(any(yap1_mamml$inframe))
})
