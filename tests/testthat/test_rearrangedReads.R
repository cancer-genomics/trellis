context("rearrangedReads")
##
## A consensus sequence for a rearrangement obtained from Delly
##
## - this unit test checks the consensus sequnce by BLAT against the hg19 reference genome
##

test_that("split reads", {
  extdata <- system.file("extdata", package="svbams")
  blat <- readBlat(file.path(extdata, "consensus.txt"))
  linked_bins <- readRDS(file.path(extdata, "consensus_linkedbins.rds"))
  lb <- linked_bins[1]
  lb$linked.to <- linked_bins[2]
  names(lb) <- "test.rid"
  expect_is(lb, "GRanges")

  ##
  ## rearrangedReads
  ##
  linked_bins <- lb
  maxgap <- 500
  is_na <- is.na(blat$Tstart)
  if(any(is_na)){
    blat <- blat[!is_na, ]
  }
  blat <- blat[blat$Tname %in% seqlevels(linked_bins), ]
  blat <- blat %>%
    mutate(match=match/Qsize*100)
  blat_gr <- blat_to_granges(blat, seqinfo(linked_bins))
  expect_identical(nrow(blat), length(blat_gr))
  genome(blat_gr) <- genome(linked_bins)
  ##
  ## A blat record must overlap one of the intervals in a linked bin
  ##
  is.overlap <- overlapsLinkedBins(blat_gr, linked_bins, maxgap=maxgap)
  blat_gr <- blat_gr[ is.overlap ]
  expect_identical(length(blat_gr), 1L)
  ## We are looking for split read alignments--a read with a blockcount of 1
  ## must have at least 2 alignments. Get rid of all reads with only a single
  ## alignment having a block count of 1. Among these reads, remove alignments that have a match score greater than 95.
  ##
  ## Alternatively, blat can report a single alignment with multiple blocks.
  ## Here, the high match score should be high, two of the blocks should hit the linked bins, and there should be more than start site in the target genome.
  ##
  tstarts <- integer_vector(blat_gr$tStarts)
  is_sr_candidate <- (blat_gr$match < 95 & blat_gr$blockcount == 1) |
    (blat_gr$match > 95 & blat_gr$blockcount > 1 & length(tstarts) > 1)
  expect_true(is_sr_candidate)
  expect_true( candidateSplitRead(blat_gr) )
  ##
  ## multipleAlignmentRecords2
  ##
  records <- blat_gr[ candidateSplitRead(blat_gr) ]
  expect_true(length(records) == 1)
})


test_that("consensus", {
  extdata <- system.file("extdata", package="svbams")
  blat <- readBlat(file.path(extdata, "consensus.txt"))
  linked_bins <- readRDS(file.path(extdata, "consensus_linkedbins.rds"))
  lb <- linked_bins[1]
  lb$linked.to <- linked_bins[2]
  names(lb) <- "test.rid"
  expect_is(lb, "GRanges")
  test <- rearrangedReads(lb, blat)
  expect_true(length(test) == 1)
##
##  blat_gr <- blat_to_granges(blat, lb)
##  expect_true("tStarts" %in% colnames(mcols(blat_gr)))
##  expect_is(blat_gr, "GRanges")
##  ##
##  ## Remove any blat alignment that does not overlap the
##  ## candidate rearrangement
##  ##
##  is.overlap <- overlapsLinkedBins(blat_gr, lb)
##  expect_is(is.overlap, "logical")
##  blat_gr <- blat_gr[is.overlap]
##
##  start.list <- blatStartList(blat_gr)
##  expect_identical(start.list,
##                   list(c(209945665L, 209947429L, 209947439L)))
##  ##
##  ## Candidate split reads either have a low match score
##  ## (only part of the read is aligned), or there are multiple
##  ## starts listed and the overall match is high
##  ##
##  is.candidate <- candidateSplitRead(blat_gr)
##  expect_false(is.candidate)
##
##  ## replace with dorothy example
##
##  if(FALSE){
##    ## list the blat alignments by tag id (qname)
##    blat_gr <- blat_gr[ is.candidate ]
##    gr.list <- split(blat_gr, blat_gr$qname)
##    n.split <- sapply(gr.list, numberAlignmentRecords)
##    expect_identical(as.integer(n.split), 3L)
##
##    gr <- multipleAlignmentRecords2(blat_gr)
##    expect_identical(length(gr), 1L)
##
##    overlaps_both <- overlapsBlatRecord(lb, gr, maxgap=500)
##    expect_true(overlaps_both)
##
##    record_overlaps <- overlapsAny(gr, lb, maxgap=500) |
##      overlapsAny(gr, lb$linked.to, maxgap=500)
##    expect_true(record_overlaps)
##    records <- gr[record_overlaps]
##    ##
##    ## Assign each sequence (qname) to a rearrangement
##    ##
##    rear.id <- get_rearrangement_id(records, lb)
##    expect_identical(rear.id, "test.rid")
##    records$rear.id <- rear.id
##    ##
##    ## split records by qname
##    ##
##    records_qname <- split(records, records$qname)
##    ##
##    ## Each sequence should have a unique rearrangement
##    ##
##    n.rear <- sapply(records_qname, function(x) length(unique(x$rear.id)))
##    records_qname <- records_qname [ n.rear == 1 ]
##    expect_identical(records_qname, GRangesList(consensus=records))
##    ##
##    ## we should check this after reduce
##    ##
##    blat.grl <- eachBlockAsGRanges(records_qname)
##    expect_true(all(elementNROWS(blat.grl) == 2))
##    blat.grl <- blat.grl[ elementNROWS(blat.grl) == 2 ]
##    ##
##    ## The split should involve nearly non-overlapping subsequences of
##    ## the read, and the total match score should be high (near 100)
##    ##
##    splitread_ranges <- lapply(blat.grl, sequenceRanges2)
##    splitread_int <- splitreadIntersection(splitread_ranges)
##    ##
##    ## total score should be near 100
##    splitread_match <- as.integer(sapply(splitread_ranges, function(x){
##      min(sum(x$match), x$Qsize)
##    }))
##    p <- splitread_match/blat_gr$Qsize[1] * 100
##    is_splitread <- splitread_int < 10 & p > 90
##    expect_true(is_splitread)
##
##    blat.grl <- blat.grl[ is_splitread ]
##    ##records_qname <- records_qname[ splitread_int < 10 & splitread_match > 90 ]
##    ##records_qname <- GRangesList(records_qname)
##    ##records <- unlist(records_qname)
##    records <- unlist(blat.grl)
##    names(records) <- NULL
##    records_rear <- split(records, records$rear.id)
##
##    # Adding seqinfo to records_rear from records_qname because rearrangedReads now
##    # maintains proper seqInfo
##    seqlevels(records_rear) <- seqlevels(records_qname)
##    seqlengths(records_rear) <- seqlengths(records_qname)
##    genome(records_rear) <- genome(records_qname)
##
##    expect_identical(length(records_rear[[1]]), 2L)
##
##    expect_identical(records_rear, test)
##  }
})

##
## A denovo deletion from paired end reads was identified on chr1
##
## - this unit test check for unmapped reads with a mapped mate near the putative sequence junction can be aligned via a split-read with blat
##
test_that("unmapped_reads_near_consensus", {
  extdata <- system.file("extdata", package="svbams")
  unmap.file <- file.path(extdata, "blat_mapped-unmapped-oralcleft.txt")
  blat <- readBlat(unmap.file)
  ##
  ## - the split reads might have been mapped by gatk
  ##
  ## - the linked bin region is not big enough for findings reads with unmapped mates
  ##
  linked_bins <- readRDS(file.path(extdata, "consensus_linkedbins.rds"))
  lb <- linked_bins[1]
  lb$linked.to <- linked_bins[2]
  expect_is(lb, "GRanges")
  names(lb) <- "test.rid"
  ## there are no split reads in this blat file
  blat_gr <- blat_to_granges(blat, seqinfo(lb))
  split_reads <- rearrangedReads2(lb, blat_gr)
  ## the split reads do not map uniquely to both ends of the target sequence
  expect_true(is.null(split_reads))
})

test_that("rearrangedReadsFun", {
  extdata <- system.file("extdata", package="svbams")
  unmap.file <- file.path(extdata, "blat_unmapped.txt")
  blat_unmap <- readBlat(unmap.file)
  rlist <- readRDS(file.path(extdata, "rlist_cgov44t.rds"))
  ##trace(rearrangedReads, browser)
  split_reads <- rearrangedReads(linkedBins(rlist), blat_unmap, 500)
  expect_identical(names(split_reads), c("1-2", "3-4"))
  qnms1 <- unique(split_reads[[1]]$qname)
  if(FALSE){
    saveRDS(qnms1, file="rearrangedReads.fbb19b6.rds")
    expected <- readRDS(file.path(extdata, "rearrangedReads.fbb19b6.rds"))
    expect_true(all(qnms1 %in% expected))
  }
  if(FALSE){
    split_reads2 <- rearrangedReadsFromRlist(rlist, blat_unmap, 500)
    expect_identical(names(split_reads2), names(split_reads))
    qnms2 <- unique(split_reads2[[1]]$qname)
    expect_identical(qnms2, qnms1)
    expect_identical(as.integer(elementNROWS(split_reads2)),
                     c(60L, 12L))
  }
})

##test_that("overlapsBoth", {
##  data(rearrangement_list)
##  data(blat_unmapped)
##  data(rearrangement_list)
##  rearranged_reads <- rearrangedReads(rearrangement_list, blat_unmapped, 500)
##  expect_is(rearranged_reads, "GRangesList")
##  ## only one rearrangement supported
##  expect_identical(length(rearranged_reads), 1L)
##  expect_true(length(rearrangement_list[names(rearranged_reads)])==1)
##  ## number split reads
##  n_splitreads <- length(split(rearranged_reads[[1]], rearranged_reads[[1]]$qname))
##  expect_identical(n_splitreads, 32L)
##
##  ## this should exit gracefully with no rearranged reads identified
##  rearranged_reads <- rearrangedReads(rearrangement_list, blat_unmapped[1:10, ], 500)
##  expect_identical(length(rearranged_reads), 0L)
##  expect_is(rearranged_reads, "GRangesList")
##})
