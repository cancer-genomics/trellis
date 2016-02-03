context("BLAT scores")
test_that("blatScores", {
  ## simulated blat record
  qnames <- c(paste0(letters[1:10], "_R1"),
              paste0(letters[1:10], "_R2"))
  ## only 1 location
  numberAlignedLocations <- rep(1, length(qnames))
  matchScores <- rep(95, length(qnames))
  Tend <-  150
  Tstart <- 51
  sblat <- data.frame(Qname=qnames,
                      match=matchScores,
                      Tend=Tend,
                      Tstart=Tstart,
                      Tname=rep("chr1", length(qnames)),
                      strand=rep("+", length(qnames)),
                      stringsAsFactors=FALSE)
  stags <- data.frame(qname=rep(letters[1:10], 2),
                      read=rep(c("R1", "R2"), each=10),
                      seqnames=rep("chr1", 10),
                      strand=rep("+", 10),
                      start=Tstart,
                      end=Tend,
                      seq=replicate(10, paste(sample(c("g","c"), 10, replace=TRUE), collapse="")),
                      stringsAsFactors=FALSE)
  stags$id <- "CGOV32T"
  stags$rearrangement.id <- "1-3"
  ## For each tag, calculate the number of near-perfect matches of the
  ## right size that overlap with eland.  If the number is 0 or more
  ## than 1, then the tag 'fails'.
  blat <- svalignments:::annotateBlatRecords(sblat, stags)
  expect_true(nrow(blat) == 20)
  s <- blatScores(sblat, stags, "SOME_ID")
  expect_true(all(s$passQC))

  ## Add reads that map to different locations with low alignment scores
  matchScores <- rep(95, length(qnames))
  Tend <-  150
  Tstart <- 51
  tmp <- data.frame(Qname=qnames,
                    match=matchScores-50,
                    Tend=Tend,
                    Tstart=Tstart,
                    Tname=rep("chr2", length(qnames)),
                    strand="+",
                    stringsAsFactors=FALSE)
  sblat2 <- rbind(sblat, tmp)
  s2 <- blatScores(sblat2, stags, "SOME_ID")
  expect_true(all(s2$passQC))

  ## Add reads that map to different locations with high alignment
  ## scores (should fail because none pass)
  tmp$match <- rep(91, length(qnames))
  sblat3 <- rbind(sblat, tmp)
  s3 <- blatScores(sblat3, stags, "SOME_ID")
  expect_true(!any(s3$passQC))

  ## Check that if 90 percent of the reads have one unique alignment
  ## that it still passes QC
  tmp$match <- rep(80, length(qnames))
  tmp$match[1] <- 90
  sblat3 <- rbind(sblat, tmp)
  s3 <- blatScores(sblat3, stags, "SOME_ID")
  expect_true(all(s3$passQC))

  ## Check that if there is one perfect match for each tag, but none
  ## overlap ELAND that it fails
  tmp$match <- rep(10, length(qnames))
  sblat4 <- rbind(sblat, tmp)
  stags4 <- stags
  stags4$seqnames=rep("chrX", 10)
  s4 <- blatScores(sblat4, stags4, "SOME_ID")
  expect_true(!any(s4$passQC))
})
