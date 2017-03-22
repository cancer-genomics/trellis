context("Subset improper read pairs")

.test_that <- function(name, expr) NULL

.test_that("addImproperReadPairs2", {
  library(svbams)
  bamdir <- system.file("extdata", package="svbams")
  bamfile <- file.path(bamdir, "cgov10t.bam")
  irpfile <- file.path(bamdir, "cgov10t_irp.rds")
  irp.all <- readRDS(irpfile)

  chroms <- c("chr3", "chr3", "chr3", "chr4", "chr9", "chrX")
  starts <- c(59600000, 60246000, 175025000, 12660000, 9800000, 6686000)
  ends <- c(61000000, 60700000, 175167000, 12770000, 10500000, 7000000)
  bins <- GRanges(chroms, IRanges(starts, ends))
  bview <- BamViews(bamPaths=bamfile, bamRanges=bins)
  aview <- AlignmentViews2(bview, irpfile)

  cbsfile <- file.path(bamdir, "cgov10t_cbs.rds")
  gr <- readRDS(cbsfile)
  gr <- gr[gr$seg.mean < log2(0.75)]
  irp <- addImproperReadPairs2(gr, aview)
  seqlevels(irp) <- seqlevels(gr)
  seqinfo(irp) <- seqinfo(gr)

  irp2 <- improperRP(gr, irp.all)
  expect_equivalent(irp, irp2)
})
