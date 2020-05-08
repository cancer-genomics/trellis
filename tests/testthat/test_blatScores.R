context("BLAT scores")
test_that("blatScores", {
  library(tidyr)
  library(tibble)
  ## simulated blat record
  qnames <- c(paste0(letters[1:10], "_R1"),
              paste0(letters[1:10], "_R2"))
  ## only 1 location
  numberAlignedLocations <- rep(1, length(qnames))
  matchScores <- rep(95, length(qnames))
  Tend <-  150
  Tstart <- 51
  ## example blat record
  sblat <- tibble(Qname=qnames,
                  match=matchScores,
                  Tend=Tend,
                  Tstart=Tstart,
                  Tname=rep("chr1", length(qnames)),
                  strand=rep("+", length(qnames)),
                  Qsize=100,
                  blockcount=1)
  ## The tags that were aligned by blat and their alignment by the whole-genome-aligner
  stags <- tibble(qname=rep(letters[1:10], 2),
                  read=rep(c("R1", "R2"), each=10),
                  seqnames=rep("chr1", 20),
                  strand=rep("+", 20),
                  start=Tstart,
                  end=Tend,
                  seq=replicate(20, paste(sample(c("g","c"), 100, replace=TRUE), collapse="")),
                  id="CGOV32T",
                  rearrangement.id="1-3")
  ## For each tag, calculate the number of near-perfect matches of the
  ## right size that overlap with the whole-genome-aligner.  
  blat <- annotateBlatRecords(sblat, stags)
  expect_true(nrow(blat) == 20)
  rear.pass <- blatScores(blat, stags, "SOME_ID")
  expect_true(nrow(rear.pass) == 1)
  expect_true(all(rear.pass$passQC))
  ##is.pass <- blatStatsPerTag(blat, tag_length=10)
  ##expect_true(!any(is.pass))
  ##is.pass <- blatStatsPerTag(blat, tag_length=100)
  ##expect_true(all(is.pass))

  ## Add reads that map to different locations with low alignment scores
  matchScores <- rep(95, length(qnames))
  Tend <-  150
  Tstart <- 51
  tmp <- tibble(Qname=qnames,
                match=matchScores-50,
                Tend=Tend,
                Tstart=Tstart,
                Tname=rep("chr2", length(qnames)),
                strand="+",
                Qsize=100,
                blockcount=1)
  sblat2 <- rbind(sblat, tmp)
  ##
  ## these should all pass since off-target have low scores
  ##
  blat <- annotateBlatRecords(sblat2, stags)
  s2 <- blatScores(blat, stags, "SOME_ID")
  expect_true(all(s2$passQC))

  ## Add reads that map to different locations with high alignment
  ## scores (should fail because none pass)
  tmp$match <- rep(91, length(qnames))
  sblat3 <- rbind(sblat, tmp)
  blat <- annotateBlatRecords(sblat3, stags)
  s3 <- blatScores(blat, stags, "SOME_ID")
  expect_true(!any(s3$passQC))

  ## Check that if 90 percent of the reads have one unique alignment
  ## that it still passes QC
  tmp$match <- rep(80, length(qnames))
  tmp$match[1] <- 90
  sblat3 <- rbind(sblat, tmp)
  blat <- annotateBlatRecords(sblat3, stags)
  s3 <- blatScores(blat, stags, "SOME_ID")
  expect_true(all(s3$passQC))

  ## Check that if there is one perfect match for each tag, but none
  ## overlap the whole-genome-aligner that it fails
  tmp$match <- rep(10, length(qnames))
  sblat4 <- rbind(sblat, tmp)
  stags4 <- stags
  stags4$seqnames=rep("chrX", 20)
  blat <- annotateBlatRecords(sblat4, stags4)
  s4 <- blatScores(blat, stags4, "SOME_ID")
 expect_true(!any(s4$passQC))
})
