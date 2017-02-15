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
  bview <- BamViews(bamRanges=bins,
                    bamPaths=file.path(path, "cgov10t.bam"))

  bins$cnt <- binCounts(bview)
  bins$std_cnt <- binNormalize(bins)
  set.seed(123)
  gc.adj <- binGCCorrect(bins)
  ##
  ##  Residuals are not centered at zero ( possibly, because this is such a
  ##  small region )
  ##
  gc.adj <- gc.adj - 0.6
  gc.adji <- as.integer(round(1000*gc.adj, 0))

  preprocessdir <- tempdir()
  fname <- file.path(preprocessdir, rdsId(bview))
  saveRDS(gc.adji, file=fname)
  pview <- PreprocessViews2(bview)
  paths(pview) <- fname
  setScale(pview) <- 1000
  ## a preprocess views object
  dat <- CNAObject(pview)
  expect_is(dat, "CNA")

  select <- 1:10
  dat2 <- CNAObject(pview[select, ])
  expect_identical(dat[select, ], dat2)

  bins$gc.adj <- gc.adj
  bins$id <- colnames(pview)[1]
  dat3 <- CNAObject(bins, "gc.adj")
  expect_equal(dat$cgov10t.bam, dat3$cgov10t.bam, tolerance=0.001)

  seg.params <- SegmentParam()
  bins$adjusted <- bins$gc.adj
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
  g.cnv <- g[g$seg.mean < -0.5]
  del.params <- DeletionParam()
  ## find all improper reads
  ## TODO:  refactor sv_deletions
  ## - do not require AlignmentViews
  ## create an AlignmentViews object
  library(svalignments)
  iparams <- improperAlignmentParams()
  file <- paste0(tempdir(), "cgov10t.rds")
  irp <- getImproperAlignmentPairs(bview, param=iparams,
                                   mapq_thr=30,
                                   use.mcols=TRUE)
  saveRDS(irp, file=file)
  aview <- AlignmentViews2(bview, file)

  ## create a PreprocessViews object
  pview <- PreprocessViews2(bview)
  gc.adji <- as.integer(round(1000*gc.adj, 0))
  file <- file.path(tempdir(), "tmp.rds")
  saveRDS(gc.adji, file=file)
  pview@bamPaths <- file
  ##
  ## critical to set scale. otherwise, the granges_copynumber is ridiculous
  ##
  ## TODO: All of these segments but the homozygous deletion gets filtered
  ## because of the fold-change context. Would be better to only compare to the
  ## segment means of the adjacent segment. The drawback is multiple segments
  ## that have about the same mean within a deletion.
  ##
  setScale(pview) <- 1000L
  egr <- expandGRanges(g.cnv, 15*width(g.cnv))
  fc_context <- granges_copynumber(g.cnv, pview)
  expect_equivalent(as.numeric(fc_context),
                    c(-1.0480, -1.1000, -8.9300, -1.0145, -1.0520))
  fc_context2 <- granges_copynumber(egr, pview)
  expect_equivalent(as.numeric(fc_context2),
                    c(-0.8965, -0.8965, -0.8965, -0.8965, -0.8965 ))
  filters <- reduceGenomeFilters(germline_filters)
  g2 <- germlineFilters(cnv=g.cnv,
                        germline_filters=filters,
                        pview=pview)
  ##
  ## Checking against the current behavior. However, 5 CNVs would be the correct
  ## result for this region.
  ##
  expected <- GRanges("chr3", IRanges(60319001, 60434001))
  tmp <- deletion_call(aview=aview, pview=pview,
                       cnv=g,
                       gr_filters=filters)
  expect_identical(calls(tmp), "homozygous")
  sv_dels <- sv_deletions(gr=g,
                          aview=aview,
                          bview=bview,
                          pview=pview,
                          gr_filters=filters,
                          param=del.params)
})
