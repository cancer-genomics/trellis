test_that("CNAObject", {
  library(svbams)
  library(svfilters.hg19)
  library(Rsamtools)
  library(svpreprocess)
  path <- system.file("extdata", package="svbams", mustWork=TRUE)
  data(bins1kb)
  data(germline_filters, package="svfilters.hg19")
  ## normalize bin counts
  bins <- keepSeqlevels(bins1kb, "chr3", pruning.mode="coarse")
  bins <- subsetByOverlaps(bins, GRanges("chr3", IRanges(59600000, 61000000)))
  bam.file <- file.path(path, "cgov10t.bam")
  bview <- BamViews(bamRanges=bins, bamPaths=bam.file)

  bins$cnt <- binCounts(bview)
  bins$std_cnt <- binNormalize(bins)
  set.seed(123)
  gc.adj <- binGCCorrect(bins)
  ##
  ##  Residuals are not centered at zero ( possibly, because this is such a
  ##  small region )
  ##
  gc.adj <- gc.adj - 0.6
  ##gc.adji <- as.integer(round(1000*gc.adj, 0))
  bins$log_ratio <- gc.adj
  seg.params <- SegmentParam()
  bins$adjusted <- bins$log_ratio
  g <- segmentBins(bins, seg.params)
  starts <- c(59599001,
              59812001,
              60141001,
              60247001)
  ends <- c(59811001,
            60140001,
            60246001,
            60318001)
  expected <- GRanges(rep("chr3", 4),
                      IRanges(starts, ends),
                      seg.mean=c(-0.0572, -1.0646, -0.0319, -1.0862))
  expect_equivalent(head(g, 4), expected)
  if(FALSE){
    library(ggplot2)
    df <- data.frame(lr=bins$adjusted,
                     start=start(bins))
    df.segs <- as.data.frame(g)
    ggplot(df, aes(start, lr)) +
      geom_point(size=0.5, col="gray") +
      geom_segment(data=df.segs, aes(x=start, xend=end,
                                     y=seg.mean, yend=seg.mean))
  }
})
