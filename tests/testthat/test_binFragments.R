test_that("binFragments", {
  library(Rsamtools)
  library(svbams)
  library(GenomicAlignments)
  data(bins)
  ### Read in bam file and great GAlignmentPairs object
  bam_path <- system.file("extdata", package="svbams")
  bamfile <- file.path(bam_path, "cgov10t.bam")
  bviews <- BamViews(bamRanges=bins, bamPaths=bamfile)
  bfile <- BamFile(bamPaths(bviews))
  galp <- readGAlignmentPairs(bfile)[1:4]
  if(FALSE){
    ## TODO: this requires some more thought. The error returned is very informative
    glist <- grglist(galp)
    gr <- GRanges(seqnames = c("chr3", "chr3"),
                  IRanges(start = c(15586441,15586441),
                          end = c(18000000, 20000000)))
    answer <- 2:3
    result <- binFragments(galp, gr)
    expect_identical(result, answer)
  }
})
