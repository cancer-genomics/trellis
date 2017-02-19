context("Germline filters for deletions")

test_that("germlineFilters", {
  library(svfilters.hg19)
  dparam <- DeletionParam()
  filters <- reduceGenomeFilters(germline_filters)

  ddir <- system.file("extdata", package="svpreprocess",
                      mustWork=TRUE)
  cov.file <- file.path(ddir, "preprocessed_coverage.rds")
  data(pview, package="svpreprocess")
  paths(pview) <- cov.file

  data(segments, package="svcnvs")
  segs <- keepSeqlevels(segments, "chr15", pruning.mode="coarse")
  ##
  ## We only have a small region of the bam file
  ## - subset the segments
  data(deletion, package="svcnvs")
  gr <- variant(deletion)
  gr <- expandGRanges(gr, 10000)
  segs <- subsetByOverlaps(segs, gr)

  cnvs <- germlineFilters(cnv=segs,
                          germline_filters=filters,
                          pview=pview,
                          param=dparam)
  if(FALSE){
    saveRDS(cnvs, file="germlineFilters.9492f3f.rds")
  }
  cnvs.9492f3f <- readRDS("germlineFilters.9492f3f.rds")
  expect_identical(cnvs, cnvs.9492f3f)

  expected <- cnvs
  ##
  ## germlineFilters
  ##
  cnv <- segs
  not_germline <- isNotGermline(cnv, filters, dparam)
  ##
  ## idea is to compute the segment means for a much larger region and then to
  ##  calculate the fold change from the original cnv to the broader region -
  ##  for this unit test, the result should be 0
  egr <- expandGRanges(cnv, 0 * width(cnv))
  expect_identical(start(egr), start(cnv))
  cn <- granges_copynumber(egr, pview)
  expect_equal(cnv$seg.mean, as.numeric(cn),
               scale=1, tolerance=0.5)

  egr <- expandGRanges(cnv, 15 * width(cnv))
  cn_context <- granges_copynumber(egr, pview)
  fc <- (cnv$seg.mean - cn_context)
  expect_equivalent(fc < -3, c(FALSE, TRUE, FALSE))

  is_big <- isLargeHemizygous(cnv, param)
  expect_true(!any(is_big))
  select <- !is_big & not_germline & fc < -0.515
  cnv <- cnv[select]
  cnvr <- reduce(cnv)
  K <- width(cnvr) > 2000
  cnvr <- cnvr[K]
  cnv <- cnv[K]
  names(cnv) <- paste0("sv", seq_along(cnv))
  expect_identical(cnv, expected)
})
