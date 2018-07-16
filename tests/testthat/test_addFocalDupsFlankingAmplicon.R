context("Adding focal duplications that flank amplicons")

test_that("addFocalDupsFlankingAmplicon", {
  ##
  ## We already have a regression test for this function in amplicon_regressions
  ##
  ## - these unit tests breakdown this function into its components
  ##
  extdata <- system.file("extdata", package="svbams")
  bamfile <- file.path(extdata, "cgov44t_revised.bam")
  path <- system.file("extdata", package="svbams")
  ag <- readRDS(file.path(path, "makeAGraphffab104.rds"))
  params <- ampliconParams()
  merged <- joinNearGRanges(ranges(ag), params)
  names(merged) <- ampliconNames(merged)
  ## the names of the nodes no longer correspond to the range names
  ## stopifnot(nodes(ag) %in% names(tmp)), and so
  ## setAmpliconGroups fails
  ranges(ag) <- merged
  rp <- get_readpairs(ag, bamfile)
  flanks <- flankingDuplications(ag, minimum_foldchange=params[["LOW_THR"]])
  ## these are two segments that are amplified to a lesser extent than the most
  ## amplified amplicon on chr8
  left.flanks <- c(flanks$left$first, flanks$left$last)
  if(FALSE){
    saveRDS(flanks, file="flankingDuplications.21adcf5.rds")
  }
  path <- system.file("extdata", package="svbams")
  flanks.21adcf5 <- readRDS(file.path(path, "flankingDuplications.21adcf5.rds"))
  flanks.21adcf5 <- updateObject(flanks.21adcf5)
  expect_identical(flanks, flanks.21adcf5)
  names(left.flanks) <- NULL
  expected <- GRanges("chr8", IRanges(128515001, 128691001),
                      seqinfo=seqinfo(merged))
  expect_identical(granges(left.flanks[1]), expected)
  expected <- GRanges("chr8", IRanges(128692001, 129163001),
                      seqinfo=seqinfo(merged))
  expect_identical(granges(left.flanks[2]), expected)
  ## nothing on the right
  right1 <- granges(flanks$right$first)
  right2 <- granges(flanks$right$last)
  if(FALSE){
    flanks.df <- as.data.frame(left.flanks)
    segs <- as.data.frame(ranges(ag))
    segs <- segs[!is.na(segs$seg.mean), ]
    ggplot(segs, aes(x=segs$start,
                     xend=segs$end,
                     y=segs$seg.mean,
                     yend=segs$seg.mean)) +
      geom_segment(aes(color=is_amplicon), size=3) +
      scale_color_manual(values=c("gray", "blue")) +
      facet_wrap(~seqnames)
    ## show the flanking regions
    segs.chr8 <- segs[segs$seqnames=="chr8", ]
    ggplot(segs.chr8, aes(x=start,
                          xend=end,
                          y=seg.mean,
                          yend=seg.mean)) +
      geom_segment(aes(color=is_amplicon), size=3) +
      geom_segment(data=flanks.df,
                   aes(x=start,
                       xend=end,
                       y=seg.mean,
                       yend=seg.mean),
                   inherit.aes=FALSE, color="black",
                   size=2)
  }
  rpsegs <- readPairsAsSegments(rp)
  ##
  ## linkedDuplicatedRanges
  ##
  flank <- flanks
  lengths <- unlist(lapply(flank, elementNROWS))
  ## this is the gap between the two left flanking ranges
  gapsLeft <- GRanges(seqnames(flank[["left"]]$first),
                      IRanges(end(flank[["left"]]$first)+1,
                              start(flank[["left"]]$last)))
  gapsRight <- GRanges(seqnames(flank[["right"]]$first),
                       IRanges(end(flank[["right"]]$first)+1,
                               start(flank[["right"]]$last)))
  cntsLeft <- minimumBasepairCoverage(rpsegs, gapsLeft)
  cntsRight <- minimumBasepairCoverage(rpsegs, gapsRight)
  minimum_count <- params[["minimum_count"]]
  left <- lapply(flank$left, "[", which(cntsLeft >= minimum_count))
  right <- lapply(flank$right, "[", which(cntsRight >= minimum_count))
  dup_gr <- list(left=left, right=right)##flanking_duplications[cnts >= minimum_count]
  ## does nothing since dup_gr is length 0
  ag2 <- addFlanks(ag, dup_gr)
  expect_identical(ag2, ag)
})
