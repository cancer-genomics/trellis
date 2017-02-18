context("Adding focal duplications that flank amplicons")

test_that("addFocalDupsFlankingAmplicon", {
  extdata <- system.file("extdata", package="svbams")
  bamfile <- file.path(extdata, "cgov44t_revised.bam")
  ag <- readRDS("makeAGraphffab104.rds")
  merged <- joinNearGRanges(ranges(ag), params)
  names(merged) <- ampliconNames(merged)
  ## the names of the nodes no longer correspond to the range names
  ## stopifnot(nodes(ag) %in% names(tmp)), and so
  ## setAmpliconGroups fails
  ranges(ag) <- merged
  rp <- svalignments::get_readpairs(ag, bamfile)

  flanks <- flankingDuplications(ag, minimum_foldchange=params[["LOW_THR"]])

  ## these are two segments that are amplified to a lesser extent than the most
  ## amplified amplicon on chr8
  left.flanks <- c(flanks$left$first, flanks$left$last)
  if(FALSE){
    
  }
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
  left.dups <- c(left1, left2)
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

  flank <- flanks
  lengths <- unlist(lapply(flank, elementNROWS))
  gapsLeft <- GRanges(seqnames(flank[["left"]]$first),
                      IRanges(end(flank[["left"]]$first)+1,
                              start(flank[["left"]]$last)))

  gapsRight <- GRanges(seqnames(flank[["right"]]$first),
                       IRanges(end(flank[["right"]]$first)+1,
                               start(flank[["right"]]$last)))
  cntsLeft <- minimumBasepairCoverage(rpsegs, gapsLeft)
  cntsRight <- minimumBasepairCoverage(rpsegs, gapsRight)
  left <- lapply(flank$left, "[", which(cntsLeft >= minimum_count))
  right <- lapply(flank$right, "[", which(cntsRight >= minimum_count))
  

  dup_gr <- linkedDuplicatedRanges(ag, rpsegs, flanks)


  L <- max(length(dup_gr$left[[1]]), length(dup_gr$right[[2]]))
  object <- addFlanks(object, dup_gr)


  ag <- addFocalDupsFlankingAmplicon(ag, rp, LOW_THR)
})
