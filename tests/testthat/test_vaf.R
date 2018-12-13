context("Variant allele frequency")

test_that("VAF for rearrangements", {
  ##
  ## Load previously identified rearrangements
  ##
  skip("Requires ovarian.cell.line.data package")
  library(svfilters.hg19)
  library(svbams)
  library(ovarian.cell.line.data)
  data(rearrangements_unordered)
  r <- rearrangements_unordered[["CGOV44T"]]
  lb <- linkedBins(r)
  r.chr8 <- r[ chromosome(lb) == "chr8" ]

  ##
  ## Annotate rearrangements
  ##




  nsplit <- elementNROWS(splitReads(r.chr8))
  nimproper <- elementNROWS(lapply(r.chr8, improper))
})
