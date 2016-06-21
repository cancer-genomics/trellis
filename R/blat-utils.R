#' Read output from the command-line blat
#'
#' The command-line version of blat can be downloaded from sourceforge:
#' \url{http://sourceforge.net/projects/blat/files/}
#' 
#' @export
#' @param filename character string providing full path to blat output
#' @return a \code{data.frame} of alignment records from blat
#' 
#' @seealso See \code{\link{blatScores}} for evaluating whether the
#'   blat alignments support a novel sequence junction
#' 
readBlat <- function(filename){
  blat <- read.delim(filename, skip=2, nrows=5, stringsAsFactors=FALSE, sep="\t", header=FALSE)
  nms <- paste0(as.character(blat[1, ]), as.character(blat[2, ]))
  nms <- gsub(" ", "", nms)
  ##nms[22:23] <- c("seq1", "seq2")
  blat <- read.delim(filename, skip=5, stringsAsFactors=FALSE, header=FALSE, sep="\t")
  colnames(blat) <- nms
  blat$Tname <- gsub(".fa", "", blat$Tname)
  blat
}

blatGRanges <- function(blat, sl=paste0("chr", c(1:22, "X", "Y", "M"))){
  g <- GRanges(blat$Tname, IRanges(blat$Tstart, blat$Tend),
               strand=factor(blat$strand, levels=c("+", "-", "*")))
  seqlevels(g, force=TRUE) <- sl
  g
}

elandGRanges <- function(blat, sl=paste0("chr", c(1:22, "X", "Y", "M"))){
  g <- GRanges(blat$eland.chr, IRanges(blat$eland.start, blat$eland.end))
  seqlevels(g, force=TRUE) <- sl
  g
}

elandAndBlatOverlap <- function(blat){
  blat$Qname <- factor(blat$Qname, levels=unique(blat$Qname))
  g.blat <- blatGRanges(blat)
  g.eland <- elandGRanges(blat)
  if(length(g.blat) != length(g.eland)){
    msg <- "Check whether unique(blat$Tname) is the same as unique(blat$eland.chr)"
    stop(msg)
  }
  strand(g.eland) <- strand(g.blat)
  same_seqname <- chromosome(g.eland) == chromosome(g.blat)
  tmp <- pintersect(g.eland[same_seqname], g.blat[same_seqname]) ##resolve.empty="start.x")
  is_overlap <- rep(FALSE, length(g.eland))
  is_overlap[same_seqname] <- width(tmp) > 0
  is_overlap
}

annotateBlatRecords <- function(blat, tag.sequences){
  tmp <- tag.sequences
  qname2 <- paste0(tmp$qname, "_", tmp$read)
  query.sequences <- setNames(tmp$seq, qname2)
  query.sequences <- query.sequences[names(query.sequences) %in% blat$Qname]
  query.sequences <- query.sequences[blat$Qname]
  ##identical(names(query.sequences), blat$Qname)
  blat$Qsequence <- query.sequences

  qnames <- factor(blat$Qname, levels=unique(blat$Qname))
  indices <- split(1:nrow(blat), qnames)
  ## T is for target (reference genome)
  ## Q is for query (tag sequence)
  sampleids <- tag.sequences$id
  ##sampleids <- rep(sampleids, elementNROWS(tag.sequences))
  eland.starts <- tag.sequences$start
  eland.ends <- tag.sequences$end
  eland.chr <- tag.sequences$seqnames
  eland.strand <- tag.sequences$strand
  ##
  ## This no longer applies
  ##
  tagnames <- paste(tag.sequences$qname, tag.sequences$read, sep="_")
  names(sampleids) <- tagnames
  names(eland.starts) <- tagnames
  names(eland.ends) <- tagnames
  names(eland.chr) <- tagnames
  names(eland.strand) <- tagnames
  ##
  ##
  ##
  ##table(blat$Qname %in% tagnames)
  qname_in_tgs <- blat$Qname %in% tagnames
  if(!all(qname_in_tgs)){
    stop("Qnames in blat records not in tag names. One of the objects must be outdated.")
    blat <- blat[blat$Qname %in% tagnames, ]
  }
  blat$id <- sampleids[blat$Qname]
  blat$eland.chr <- eland.chr[blat$Qname]
  blat$eland.start <- eland.starts[blat$Qname]
  blat$eland.end <- eland.ends[blat$Qname]
  blat$eland.strand <- eland.strand[blat$Qname]
  ## replace Tsize with T.end-T.start
  blat$Tsize <- blat$Tend-blat$Tstart
  blat$is_overlap <- elandAndBlatOverlap(blat)
  rownames(blat) <- NULL
  blat
}

## A read can not have 2 hits meeting the 90-90 rule
removeReadsWithoutMate <- function(blat){
  rid <- blat$Qname
  rpid <- gsub("_R[12]$", "", rid)
  rpid <- factor(rpid, levels=unique(rpid))
  rplist <- lapply(split(rid, rpid), function(x) unique(x))
  el <- elementNROWS(rplist)
  ## there are 36 reads with no mate listed
  rps.without.mate <- names(rplist)[el != 2]
  blat <- blat[!rpid %in% rps.without.mate, ]
}

.listTagsByGroup <- function(tags, tag.group){
  tag.group <- factor(tag.group, levels=unique(tag.group))
  tagnames <- paste(tags$qname, tags$read, sep="_")
  tag.list <- split(tagnames, tag.group)
}

addXCoordinateForTag <- function(blat){
  qnm <- blat$Qname
  rid <- gsub("_R[12]", "", qnm)
  x <- as.integer(factor(rid, levels=unique(rid)))
  is_r2 <- rep(FALSE, length(x))
  is_r2[grep("_R2$", qnm)] <- TRUE
  ## for plotting
  x[!is_r2] <- x[!is_r2] - 0.2
  x[is_r2] <- x[is_r2] + 0.2
  x
}

blatStatsPerTag <- function(blat.records, tag_length){
  ## summary statistics for a single read
  is_90 <- blat.records$match > 90
  ## Tsize in blat.records is size of target sequence (chromosome)
  Tsize <- abs(blat.records$Tend-blat.records$Tstart)
  is_size_near100 <- Tsize > (tag_length - 1/5*tag_length) &
    Tsize < (tag_length + 1/5*tag_length)
  is_90 <- is_90 & is_size_near100
  is_overlap <- blat.records$is_overlap
  ## calculate number of near-perfect matches for each tag
  n.matches <- sapply(split(is_90, blat.records$Qname), sum)
  n.eland.matches <- sapply(split(is_overlap & is_90, blat.records$Qname), sum)
  ## there should only be one eland alignment with high quality (n.matches == 1)
  ## AND this read should overlap the eland alignment
  n.matches==1 & n.eland.matches==1
}

.blatStatsRearrangement <- function(blat, thr=0.8, tag_length){
  cols <- c("Qname", "match", "is_overlap", "Tstart", "Tend")
  blat <- blat[, cols]
  if(nrow(blat) == 0) return(NULL)
  blat$tag_index <- addXCoordinateForTag(blat)
  stats <- blatStatsPerTag(blat, tag_length)
  proportion_pass <- mean(stats)
  qnames <- gsub("R[12]", "", blat$Qname)
  n_tags <- length(unique(qnames))
  is_pass <- proportion_pass > thr & n_tags >= 5
  ##blat$passQC <- rep(proportion_pass > thr, nrow(blat))
  blat$passQC <- is_pass
  blat
}

#' blatScores assesses whether the improper read pairs at a
#' rearrangement junction provide strong support of the rearrangement
#'
#' In the following, we refer to a 'record' as one row in the table of
#' blat output -- i.e., one (of possibly many) alignments for a read.
#' This function adds an indicator for whether the blat alignment is
#' consistent with the original alignment (\code{is_overlap}), an
#' indicator of whether it passes quality control (\code{passQC}) (see
#' details), the id of the rearrangement (\code{rearrangement}), and
#' the sample id (\code{id}).
#'
#' @details
#'
#' \strong{Read-level QC:} 
#'
#' For each record, we evaluate
#'
#' \enumerate{
#' 
#' \item whether the matching score is at least 90
#' 
#' \item whether the size of the target alignment (Tstart - Tend) is
#'   less than 120bp and more than 80bp
#'
#' \item whether the blat alignment overlaps with any of the original
#'   alignment
#'
#' }
#'
#' The records are then grouped by Qname.  For each Qname, we compute
#'   the number of reads with a score above 90 and within the
#'   specified size range.  In addition, we compute the number of
#'   reads with a score above 90, within the specified size range, and
#'   that overlap with the original alignment.  The Qname passes QC if
#'   both sums evaluate to 1.  That is, a read passes QC only if there
#'   is a single record with a high BLAT score within the specified
#'   size range and this single record overlaps with the original
#'   alignment.
#'
#' \strong{Rearrangement-level QC:}
#'
#' Given a pass / fail designation for each read by the above
#'   analysis, we group the reads by the rearrangement id.  A
#'   rearrangement passes QC if > 80% of the reads supporting the
#'   rearrangement pass QC.
#' 
#' @examples
#'  qnames <- c(paste0(letters[1:10], "_R1"),
#'              paste0(letters[1:10], "_R2"))
#'  ## only 1 location
#'  numberAlignedLocations <- rep(1, length(qnames))
#'  matchScores <- rep(95, length(qnames))
#'  Tend <-  150
#'  Tstart <- 51
#'  ## Made up output from blat
#'  sblat <- data.frame(Qname=qnames,
#'                      match=matchScores,
#'                      Tend=Tend,
#'                      Tstart=Tstart,
#'                      Tname=rep("chr1", length(qnames)),
#'                      strand=rep("+", length(qnames)),
#'                      stringsAsFactors=FALSE)
#' 
#' ## Made up information on a rearrangement, including the locations
#' ##  at which the reads were originally aligned
#' 
#'  stags <- data.frame(qname=rep(letters[1:10], 2),
#'                      read=rep(c("R1", "R2"), each=10),
#'                      seqnames=rep("chr1", 10),
#'                      strand=rep("+", 10),
#'                      start=Tstart,
#'                      end=Tend,
#'                      seq=replicate(10, paste(sample(c("g","c"), 10, replace=TRUE),
#'                          collapse="")),
#'                      stringsAsFactors=FALSE)
#'  stags$id <- "CGOV32T"
#'  stags$rearrangement.id <- "1-3"
#'  ## For each tag, calculate the number of near-perfect matches of the
#'  ## right size that overlap with eland.  If the number is 0 or more
#'  ## than 1, then the tag 'fails'.
#'  blat <- svalignments:::annotateBlatRecords(sblat, stags)
#'  nrow(blat) == 20
#'  s <- blatScores(sblat, stags, "SOME_ID")
#'  head(s)
#' @export
#' @param blat a \code{data.frame} of results from  command-line blat
#' @param tags a \code{data.frame} containing read names and the original alignment locations
#' @param id  a character-vector of sample identifiers
#' @param thr a length-one numeric vector indicating the fraction of
#'   reads at a rearrangement that must pass the read-level QC.
blatScores <- function(blat, tags, id, thr=0.8){
  tag_length <- nchar(tags$seq[1])
  blat$match <- blat$match/tag_length * 100
  rownames(tags) <- paste0(tags$qname, "_", tags$read)
  blat <- annotateBlatRecords(blat, tags)
  blat <- removeReadsWithoutMate(blat)
  tagnames.list <- .listTagsByGroup(tags, tags[["rearrangement.id"]])
  blat.parsed <- vector("list", length(tagnames.list))
  for(j in seq_along(tagnames.list)){
    tagnames <- tagnames.list[[j]]
    rid <- names(tagnames.list)[j]
    blat_rid <- blat[blat$Qname %in% tagnames, ]
    result <- .blatStatsRearrangement(blat_rid, thr=thr, tag_length)
    if(!is.null(result)){
      result$rearrangement <- rep(rid, nrow(result))
    }
    blat.parsed[[j]] <- result
  }
  blat.parsed <- blat.parsed[!sapply(blat.parsed, is.null)]
  blat.parsed <- do.call("rbind", blat.parsed)
  blat.parsed$id <- rep(id, nrow(blat.parsed))
  blat.parsed
}


#' Reads files output by the blat executable.
#'
#' The blat executable is run with the default set of arguments using
#' a reference genome comprised of standard sequences without
#' alternates.  This function is a wrapper for \code{blatScores}. If a
#' file containing \code{blatScores} already exists, this function
#' only reads previously saved computations from disk.
#'
#' @seealso See \code{blatScores}
#' @examples
#'  library(svovarian)
#'  dirs <- projectOvarian(rootname="OvarianData2")
#'  if(FALSE){
#'    tags <- readRDS(file.path(dirs[["3read"]], "CGOV2T.rds"))
#'    blat <- readBlat(file.path(dirs["1blat"], "CGOV2T.txt"))
#'    parsed <- scoreBlatExperiment(id, blat, tags, dirs)
#'    parsed
#'  }
#'
#' @return a list of \code{tbl_df} objects
#' 
#' @param id sample id
#'
#' @param blat data.frame of blat records from the command-line blat tool
#' 
#' @param tags a \code{tbl_df} object of the read sequences
#' 
#' @param dirs a \code{DataPaths} object
#' 
#' @param thr a length-one numeric vector indicating the fraction of
#'   improper reads that are of high quality for a given
#'   rearrangement.
#' 
scoreBlatExperiment <- function(id, blat, tags, dirs, thr=0.8){
  parsed_file <- file.path(dirs[["4parsed_mapped"]], paste0(id, ".rds"))
  if(file.exists(parsed_file)){
    blat2 <- readRDS(parsed_file)
  } else {
    blat2 <- blatScores(blat, tags, id=id, thr=thr)
    saveRDS(blat2, file=parsed_file)
  }
  blat2 <- tbl_df(blat2)
  blat2
}

overlapsBlatRecord <- function(linked_bins, blat_record, maxgap=200){
  overlapsAny(linked_bins, blat_record, maxgap=maxgap) &
    overlapsAny(linked_bins$linked.to, blat_record, maxgap=maxgap)
}

sequenceRanges <- function(blat){
  ##
  ## We ignore the seqname, but we use GRanges instead of IRanges to
  ## make use of the metadata (here, the score from blat)
  ##
  GRanges("seq", IRanges(blat$qstart, blat$qend), match=blat$match)
}

multipleAlignmentRecords <- function(records){
  records <- records[ records$match < 95]
  recordlist <- split(records, records$qname)
  n.alignments <- elementNROWS(recordlist)
  recordlist <- recordlist[n.alignments >= 2]
  unlist(recordlist)
}

#' Identify rearranged reads -- initiallly unmapped reads that can be
#' aligned by blat to span a novel sequence junction.
#'
#' @return a list of blat records that map to both sides of a sequence
#'   junction.  Each list element corresponds to one read that is
#'   aligned to two locations.  
#' 
#' @export
#' 
#' @param r a \code{RearrangementList} object
#' 
#' @param blat a data.frame of blat alignment records
#' 
#' @param maxgap this maximum gap between the mapped read and the
#'   genomic intervals of the improper read clusters
rearrangedReads <- function(r, blat, maxgap=500){
  is_na <- is.na(blat$Tstart)
  if(any(is_na)){
    blat <- blat[!is_na, ]
  }
  lb <- linkedBins(r)
  blat <- blat[blat$Tname %in% seqlevels(lb), ]
  blat_gr <- GRanges(blat$Tname, IRanges(blat$Tstart, blat$Tend),
                     match=blat$match, qname=blat$Qname,
                     qstart=blat$Qstart, qend=blat$Qend,
                     seqinfo=seqinfo(lb))
  ##
  ## A blat record must overlap one of the intervals in a linked bin
  ##
  record_overlaps <- overlapsAny(blat_gr, lb, maxgap=maxgap) |
    overlapsAny(blat_gr, lb$linked.to, maxgap=maxgap)
  records <- blat_gr[record_overlaps]
  ## We are looking for split read alignments--a read must have at
  ## least 2 alignments.  Get rid of all reads with only a single
  ## alignment, and all alignments that have a match score greater
  ## than 95.
  records <- multipleAlignmentRecords(records)
  ##
  ## Both intervals in a linked bin must overlap a blat record
  ##
  overlaps_both <- overlapsBlatRecord(lb, records, maxgap)
  lb <- lb[overlaps_both]
  ## Repeat filter on the records
  record_overlaps <- overlapsAny(records, lb, maxgap=maxgap) |
    overlapsAny(records, lb$linked.to, maxgap=maxgap)
  records <- records[record_overlaps]
  if(length(records) == 0){
    return(GRangesList())
  }
  ##
  ## Assign each sequence (qname) to a rearrangement
  ##
  records$rear.id <- NA
  h1 <- findOverlaps(records, lb, maxgap=maxgap)
  records$rear.id[queryHits(h1)] <- names(lb)[subjectHits(h1)]
  h2 <- findOverlaps(records, lb$linked.to, maxgap=maxgap)
  records$rear.id[queryHits(h2)] <- names(lb)[subjectHits(h2)]
  ##
  ## split records by qname
  ##
  records_qname <- split(records, records$qname)  
  ##
  ## Each sequence should have a unique rearrangement
  ##
  n.rear <- sapply(records_qname, function(x) length(unique(x$rear.id)))
  records_qname <- records_qname [ n.rear == 1 ]
  ##
  ## And the number of records for a given rearrangement should be 2
  ##
  records_qname <- records_qname [ elementNROWS(records_qname) == 2 ]
  ##
  ## The split should involve nearly non-overlapping subsequences of
  ## the read, and the total match score should be high (near 100)
  ##
  splitread_ranges <- lapply(records_qname, sequenceRanges)
  ## intersection
  splitread_int <- as.integer(sapply(splitread_ranges, function(x) {
    w <- width(intersect(x[1], x[2]))
    ## return 0 if no overlap
    if(length(w) == 0) w <- 0L 
    w
  }))
  ## total score should be near 100
  splitread_match <- as.integer(sapply(splitread_ranges, function(x) sum(x$match)))
  is_splitread <- splitread_int < 10 & splitread_match > 90
  records_qname <- records_qname[ splitread_int < 10 & splitread_match > 90 ]
  records <- unlist(records_qname)
  records$qname <- names(records)
  names(records) <- NULL
  ## list the split reads by rearrangement
  records_rear <- split(records, records$rear.id)
}



#' Identify rearranged reads -- initiallly unmapped reads that can be
#' aligned by blat to span a novel sequence junction.
#'
#' 
#' For unmapped-mapped read pairs where the mapped read is aligned
#' near a candidate junction, assess whether the unmapped read spans a
#' novel sequence junctions.
#'
#' @return a list of rearranged reads
#' @export
#' @param dirs a \code{DataPaths} object
#' @param rlist a list of \code{RearrangementList} objects
#' @param maxgap a length-one integer vector specifying how much gap
#'   is allowed between a read alignment location and the genomic
#'   intervals for the rearrangement clusters
#' 
rearrangedReadList <- function(dirs, rlist, maxgap=500){
  ids <- names(rlist)
  outfiles <- file.path(dirs[["5parsed_unmapped"]], paste0(ids, ".rds"))
  infiles <- file.path(dirs[["2blat_unmapped"]], paste0(ids, ".rds"))
  results <- setNames(vector("list", length(ids)), ids)
  for(j in seq_along(ids)){
    if(file.exists(outfiles[j])){
      results[[j]] <- readRDS(outfiles[j])
    }
    blat <- readBlat(infiles[j])
    results[[j]] <- rearrangedReads(rlist[[j]], blat, maxgap)
    saveRDS(results[[j]], file=outfiles[j])
  }
  results
}


removeAmbiguous <- function(x, blat){
  blat <- blat[blat$passQC, ]
  x <- x[names(x) %in% blat$rearrangement]
  x
}

#' Removes rearrangements for which the BLAT-aligned reads do not pass QC
#'
#' BLAT records for which the \code{passQC} variable is \code{FALSE}
#' are removed.  The \code{RearrangementList} is then subset to
#' include only those rearrangement ids that remain in the BLAT
#' \code{data.frame}.
#'
#' 
#' @seealso See \code{blatScores} for the approach used to QC
#'   BLAT-aligned reads.
#' 
#' @export
#' @param rlist a \code{RearrangementList}
#' @param blat_list a list of \code{data.frame} objects, as returned
#'   by \code{scoreBlatExperiment}
#' 
#' @examples
#'  library(svovarian)
#'  dirs <- projectOvarian(rootname="OvarianData2")
#'  if(FALSE){
#'  ## A BLAT-filtered RearrangementList
#'  ##saved_result <- readRDS(file.path(dirs[["unit_test"]], "blat_filtered.rds"))
#'  ##tag_list <- readRDS(file.path(dirs[["unit_test"]], "tag_seqs.rds"))
#'  blat <- scoreBlatExperiment(tag_list, dirs)
#'  blat_list <- blat["CGOV2T"]
#'
#'  rlist <- readRDS(file.path(dirs[["rear:filter"]], "CGOV2T.rds"))
#'  rlist <- list(CGOV2T=rlist)
#'  ## Another BLAT-filtered RearrangementList created by 'removeAmbigousAln'
#'  filtered_rlist <- removeAmbiguousAln(rlist, blat_list)
#'  print(filtered_rlist)
#'  }
removeAmbiguousAln <- function(rlist, blat_list){
  mapply(removeAmbiguous, x=rlist, blat=blat_list)
}
