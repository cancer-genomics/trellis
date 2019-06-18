context("Variant allele frequency")

test_that("VAF for rearrangements", {
  ##
  ## Load previously identified rearrangements
  ##
  skip("Requires ovarian.cell.line.data package")
  library(svfilters.hg19)
  library(svbams)
  library(ovarian.cell.line.data)
  library(trellis)
  data(rearrangements_unordered)
  r <- rearrangements_unordered[["CGOV44T"]]
  lb <- linkedBins(r)
  r.chr8 <- r[ chromosome(lb) == "chr8" ]

  ##
  ## Annotate rearrangements
  ##
  library(GenomicAlignments)
  rear1 <- r.chr8[[1]]
  irp <- improper(rear1)
  r1 <- first(irp)
  r2 <- last(irp)
  reduce(as(r1, "GRanges"))
  reduce(as(r2, "GRanges"))
  ref_3prime <-   reduce(as(r1, "GRanges"))[1]
  start(ref_3prime) <- start(ref_3prime) - 300
  end(ref_3prime) <- start(ref_3prime)+300

  ref_5prime <-   reduce(as(r1, "GRanges"))[2]
  end(ref_5prime) <- end(ref_5prime) + 300
  start(ref_5prime) <- end(ref_5prime) - 300
  strand(ref_3prime) <- "*"
  strand(ref_5prime) <- "*"

  library(Rsamtools)
  scanBam(bamfile, what=c(ref_5prime, ref_3prime))

  modalRearrangement(r.chr8)



  irp <- improper(r1)
  read1 <- first(irp)
  read2 <- last(irp)

  read1.gr <- as(read1, "GRanges")
  read2.gr <- as(read2, "GRanges")
  ##
  ## Reduced representation
  ##
  ## track the index of reads that go into the reduced representation
  read1.reduced <- reduce(read1.gr)
  read2.reduced <- reduce(read2.gr)

  sr <- splitReads(r1)
  ## TODO:  Strand not recorded for split read !
  sr.reduced <- reduce(sr)
  read1.reduced
  read2.reduced
  sr.reduced





  nsplit <- elementNROWS(splitReads(r.chr8))
  nimproper <- elementNROWS(lapply(r.chr8, improper))
})
