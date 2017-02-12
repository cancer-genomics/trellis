test_that("CNAObject", {
  library(svbams)
  library(svfilters.hg19)
  library(Rsamtools)
  library(svpreprocess)
  path <- system.file("extdata", package="svbams", mustWork=TRUE)
  data(bins1kb)
  ## normalize bin counts
  bins <- keepSeqlevels(bins1kb, "chr3", pruning.mode="coarse")
  bins <- subsetByOverlaps(bins, GRanges("chr3", IRanges(59600000, 61000000)))
  bview <- BamViews(bamRanges=bins,
                    bamPaths=file.path(path, "cgov10t.bam"))

  bins$cnt <- binCounts(bview)
  bins$std_cnt <- binNormalize(bins)
  set.seed(123)
  gc.adj <- binGCCorrect(bins)
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
})
