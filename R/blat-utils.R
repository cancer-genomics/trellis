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
  seqlevels(g, pruning.mode="coarse") <- sl
  g
}

genomeGRanges <- function(blat, sl=paste0("chr", c(1:22, "X", "Y", "M"))){
  g <- GRanges(blat$genome.chr, IRanges(blat$genome.start, blat$genome.end))
  seqlevels(g, pruning.mode="coarse") <- sl
  g
}

genomeAndBlatOverlap <- function(blat){
  blat$Qname <- factor(blat$Qname, levels=unique(blat$Qname))
  g.blat <- blatGRanges(blat)
  g.genome <- genomeGRanges(blat)
  if(length(g.blat) != length(g.genome)){
    msg <- "Check whether unique(blat$Tname) is the same as unique(blat$genome.chr)"
    stop(msg)
  }
  strand(g.genome) <- strand(g.blat)
  same_seqname <- chromosome(g.genome) == chromosome(g.blat)
  tmp <- pintersect(g.genome[same_seqname], g.blat[same_seqname]) ##resolve.empty="start.x")
  is_overlap <- rep(FALSE, length(g.genome))
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

  ##
  ## genome aligner
  ##
  genome.starts <- tag.sequences$start
  genome.ends <- tag.sequences$end
  genome.chr <- tag.sequences$seqnames
  genome.strand <- tag.sequences$strand
  ##
  ## This no longer applies
  ##
  tagnames <- paste(tag.sequences$qname, tag.sequences$read, sep="_")
  names(sampleids) <- tagnames
  names(genome.starts) <- tagnames
  names(genome.ends) <- tagnames
  names(genome.chr) <- tagnames
  names(genome.strand) <- tagnames
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
  blat$genome.chr <- genome.chr[blat$Qname]
  blat$genome.start <- genome.starts[blat$Qname]
  blat$genome.end <- genome.ends[blat$Qname]
  blat$genome.strand <- genome.strand[blat$Qname]
  ## replace Tsize with T.end-T.start
  blat$Tsize <- blat$Tend-blat$Tstart
  blat$is_overlap <- genomeAndBlatOverlap(blat)
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
  stats <- blat.records %>%
    mutate(is_90 = match > 90,
           Tsize=abs(Tend-Tstart),
           is_size_near100=Tsize > (tag_length-1/5*tag_length) &
             Tsize < (tag_length + 1/5*tag_length),
           is_90 = is_90 & is_size_near100) %>%
    group_by(Qname) %>%
    summarize(n.matches=sum(is_90),
              n.eland.matches=sum(is_90 & is_overlap)) %>%
    mutate(is_pass=n.matches == 1 & n.eland.matches==1)
  return(stats)
  if(FALSE){
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
}

.blatStatsRearrangement <- function(blat, thr=0.8, tag_length){
  cols <- c("Qname", "match", "is_overlap", "Tstart", "Tend")
  blat <- blat[, cols]
  if(nrow(blat) == 0) return(NULL)
  splitby <- blat$Qname %>%
    strsplit("R[12]") %>%
    sapply("[", 1)
  blat$tag_index <- blat$Qname %>% strsplit(splitby) %>%
    sapply("[", 2)
  ##blat$tag_index <- addXCoordinateForTag(blat)
  ##stats <- blatStatsPerTag(blat, tag_length)
  stats <- blat %>%
    mutate(is_90 = match > 90,
           Tsize=abs(Tend-Tstart),
           is_size_near100=Tsize > (tag_length-1/5*tag_length) &
             Tsize < (tag_length + 1/5*tag_length),
           is_90 = is_90 & is_size_near100) %>%
    group_by(Qname) %>%
    summarize(n.matches=sum(is_90),
              n.eland.matches=sum(is_90 & is_overlap)) %>%
    mutate(is_pass=n.matches == 1 & n.eland.matches==1)

  proportion_pass <- mean(stats)
  qnames <- gsub("R[12]", "", blat$Qname)
  n_tags <- length(unique(qnames))
  is_pass <- proportion_pass > thr & n_tags >= 5
  ##blat$passQC <- rep(proportion_pass > thr, nrow(blat))
  blat$passQC <- is_pass
  blat
}

blatGRanges <- function(blat){
  chr <- blat$Tname
  starts <- blat$Tstart
  ends <- blat$Tend
  g <- GRanges(chr, IRanges(starts, ends),
               blockcount=blat$blockcount,
               Qname=blat$Qname,
               rid=blat$rid,
               score=blat$score)
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
#'  blat <- trellis:::annotateBlatRecords(sblat, stags)
#'  nrow(blat) == 20
#'  s <- blatScores(sblat, stags, "SOME_ID")
#'  head(s)
#' @export
#' @param blat a \code{data.frame} of results from  command-line blat
#' @param tags a \code{data.frame} containing read names and the original alignment locations
#' @param prop.pass a length-one numeric vector indicating the fraction of
#'   reads at a rearrangement that must pass the read-level QC.
#' @param min.tags the minimum number of tags that pass BLAT QC for each rearrangement
#' @param id sample id
blatScores <- function(blat, tags, id, min.tags=5, prop.pass=0.8){
  tags <- tags %>%
    unite("Qname", c("qname", "read"))
  ##table(blat_aln$Qname %in% tags$Qname)
  rid <- tags %>%
    group_by(Qname) %>%
    summarize(rid=unique(rearrangement.id))
  blat2 <- left_join(blat, rid, by="Qname") %>%
    mutate(score=match/Qsize) %>%
    mutate(tsize=abs(Tstart-Tend)) %>%
    filter(score >= 0.90 & blockcount==1) %>%
    filter(tsize >= (Qsize - 1/5*Qsize) & tsize <= (Qsize + 1/5*Qsize))
  blat.g <- blatGRanges(blat2)
  genome.g <- GRanges(tags$seqnames, IRanges(tags$start, tags$end),
                      Qname=tags$Qname, rid=tags$rearrangement.id)
  hits <- findOverlaps(blat.g, genome.g, ignore.strand=TRUE, maxgap=500)
  keep <- blat.g$Qname[queryHits(hits)] == genome.g$Qname[subjectHits(hits)]
  hits <- hits[keep]
  blat2$overlaps_genome <- FALSE
  blat2$overlaps_genome[ queryHits(hits) ] <- TRUE
  blat3 <- blat2 %>%
    group_by(Qname)  %>%
    summarize(overlaps_genome=sum(overlaps_genome),
              overlaps_other=sum(!overlaps_genome),
              rid=unique(rid)) %>%
    ungroup %>%
    group_by(rid) %>%
    summarize(number_reads=length(unique(Qname)),
              p_overlap_genome=mean(overlaps_genome >= 1),
              p_notoverlap_genome=mean(overlaps_genome < 1)) %>%
    mutate(passQC=p_overlap_genome >= prop.pass &
             number_reads >= min.tags) %>%
    mutate(rid=factor(rid))
  return(blat3)
  ## bind blat results with tag information
##  tag_length <- nchar(tags$seq[1])
##  blat <- blat %>%
##    mutate(score=match/tag_length*100) %>%
##    filter(score >= 90)
##
##  tags2 <- tags
##  ix <- match("qname", colnames(tags2))
##  colnames(tags2) <- paste0("genome_", colnames(tags2))
##  colnames(tags2)[ix] <- "Qname"
##  splitby <- blat$Qname %>%
##    strsplit("R[12]") %>%
##    sapply("[", 1)
##  blat$read <-   blat$Qname %>% strsplit(splitby) %>%
##    sapply("[", 2)
##  blat$Qname <- gsub("_R[12]", "", blat$Qname)
##  blat2 <- left_join(blat, tags2, by="Qname") %>%
##    as.tibble
##  blat.gr <- GRanges(blat2$Tname, IRanges(blat2$Tstart, blat2$Tend))
##  tags.gr <- GRanges(as.character(blat2$genome_seqnames),
##                     IRanges(blat2$genome_start, blat2$genome_end))
##  hits <- findOverlaps(blat.gr, tags.gr, ignore.strand=TRUE)
##  keep <- queryHits(hits) == subjectHits(hits)
##  hits <- hits[keep]
##  blat.gr$overlaps_genome <- FALSE
##  blat.gr$overlaps_genome[queryHits(hits)] <- TRUE
##  blat2$overlaps_genome <- blat.gr$overlaps_genome
##  ## Remove reads without mate
##  tab <- table(tags2$Qname)
##  tab <- tab[tab < 2]
##  if(length(tab) > 0){
##    blat2 <- filter(blat2, !Qname %in% names(tab))
##  }
##  ##blat$match <- blat$match/tag_length * 100
##  ##rownames(tags) <- paste0(tags$qname, "_", tags$read)
##  ##blat3 <- annotateBlatRecords(blat2, tags) %>%
##  ##as.tibble
##  ##blat <- removeReadsWithoutMate(blat)
##  ##colnames(tags) <- gsub("rearrangement.id", "rid", colnames(tags))
##  ix <- match("genome_rearrangement.id", colnames(blat2))
##  colnames(blat2)[ix] <- "rid"
  ##blat2 <- left_join(blat, tags, by="Qname")
##  blat2 <- .listTagsByGroup(tags, tags[["rearrangement.id"]]) %>% {
##    tibble(rid=rep(names(.), elementNROWS(.)),
##           Qname=unlist(.))
##  } %>% left_join(blat)
##  splitby <- blat2$Qname %>%
##    strsplit("R[12]") %>%
##    sapply("[", 1)
##  blat2$read <-   blat2$Qname %>% strsplit(splitby) %>%
##    sapply("[", 2)
  blat.qname <- blat2 %>%
    mutate(Tsize=abs(Tend-Tstart),
           is_size_near100=Tsize > (tag_length-1/5*tag_length) &
             Tsize < (tag_length + 1/5*tag_length),
           is_90 = match > 90 & is_size_near100) %>%
    filter(is_90) %>%
    group_by(Qname) %>%
    summarize(rid=unique(rid),
              n.matches=sum(is_90),
              n.genome.matches=sum(is_90 & overlaps_genome)) %>%
    select(c("rid", "n.matches", "n.genome.matches", "Qname")) %>%
    mutate(is_pass=n.matches == 2 & n.genome.matches  >= 2)
  blat.rid <- blat.qname %>%
    group_by(rid) %>%
    summarize(proportion_pass=mean(is_pass),
              ntags=n()) %>%
    mutate(passQC=proportion_pass > prop.pass & ntags >= min.tags)
  return(blat.rid)
  ##  blat.parsed <- vector("list", length(tagnames.list))
  ##  for(j in seq_along(tagnames.list)){
  ##    tagnames <- tagnames.list[[j]]
  ##    rid <- names(tagnames.list)[j]

  ##    blat_rid <- blat[blat$Qname %in% tagnames, ]
  ##    result <- .blatStatsRearrangement(blat_rid, thr=thr, tag_length)
  ##    if(!is.null(result)){
  ##      result$rearrangement <- rep(rid, nrow(result))
  ##    }
  ##    blat.parsed[[j]] <- result
  ##  }
  ##  blat.parsed <- blat.parsed[!sapply(blat.parsed, is.null)]
  ##  blat.parsed <- do.call("rbind", blat.parsed)
  ##  blat.parsed$id <- rep(id, nrow(blat.parsed))
  ##blat.parsed
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

sequenceRanges2 <- function(blat){
  ##
  ## We ignore the seqname, but we use GRanges instead of IRanges to
  ## make use of the metadata (here, the score from blat)
  ##
  GRanges("seq", IRanges(blat$qStarts, width=blat$blockSizes),
          match=blat$match,
          Qsize=blat$Qsize)
}

multipleAlignmentRecords <- function(records){
  records <- records[ records$match < 95]
  recordlist <- split(records, records$qname)
  n.alignments <- elementNROWS(recordlist)
  recordlist <- recordlist[n.alignments >= 2]
  unlist(GRangesList(recordlist))
}

integer_vector <- function(x){
  as.integer(unlist(strsplit(x, ",")))
}

##start_vector <- function(g$tStarts){
##  as.integer(unlist(strsplit(tStarts, ",")))
##}
##
##block_vector <- function(blockSizes){
##  as.integer(unlist(strsplit(blockSizes, ",")))
##}

blatStartList <- function(blat.gr){
  lapply(blat.gr$tStarts, integer_vector)
}

blatBlockList <- function(blat.gr){
  lapply(blat.gr$blocksizes, integer_vector)
}

candidateSplitRead <- function(blat.gr){
  blat.gr$match < 95 | blat.gr$blockcount > 1
}

numberAlignmentRecords <- function(blat.gr){
  L <- length(blat.gr)
  if(L == 1){
    L <- blat.gr$blockcount
  }
  L
}

.each_block_granges <- function(g){
  ##browser()
  starts <- integer_vector(g$tStarts)
  L <- length(starts)
  widths <- integer_vector(g$blockSizes)
  qstarts <- integer_vector(g$qStarts)
  bsizes <- integer_vector(g$blockSizes)
  chrom <- rep(chromosome(g), g$blockcount)
  qends <- qstarts+bsizes
  bmatch <- rep(g$match, g$blockcount)
  gapbases <- rep(g$gapbases, g$blockcount)
  strands <- rep(cstrand(g), g$blockcount)
  g2 <- GRanges(chrom,
                IRanges(starts, width=widths),
                qend=qends,
                qStarts=qstarts,
                blockSizes=bsizes,
                gapbases=gapbases,
                match=bmatch,
                strand=strands)
  gr <- reduce(g2, with.revmap=TRUE)
  revmap <- mcols(gr)$revmap
  tmp <- relist(g2$qStarts[unlist(revmap)], revmap)
  qstarts <- sapply(tmp, min)
  tmp <- relist(g2$blockSizes[unlist(revmap)], revmap)
  bsizes <- sapply(tmp, sum)
  tmp <- relist(g2$match[unlist(revmap)], revmap)
  bmatch <- sapply(tmp, mean)
  tmp <- relist(g2$gapbases[unlist(revmap)], revmap)
  gapbases <- sapply(tmp, min)
  L <- length(gr)
  gr$qStarts <- qstarts
  gr$blockSizes <- bsizes
  gr$match <- bmatch
  gr$rear.id <- rep(g$rear.id[1], L)
  gr$qname <- rep(g$qname[1], L)
  gr$gapbases <- gapbases
  gr$Qsize <- rep(g$Qsize[1], L)
  ##gr$qend <- ##rep(g$qend[1], L)
  gr$qend <- gr$qStarts + gr$blockSizes
  gr
}

eachBlockAsGRanges <- function(blat.grl){
##  browser()
##  for(i in seq_along(blat.grl)){
##    .each_block_granges(blat.grl[[i]])
##  }
  blat.grl2 <- lapply(blat.grl, .each_block_granges)
  GRangesList(blat.grl2)
}

multipleAlignmentRecords2 <- function(records){
  records <- records[ candidateSplitRead(records) ]
  recordlist <- split(records, records$qname)
  ##n.alignments <- elementNROWS(recordlist)
  n.alignments <- sapply(recordlist, numberAlignmentRecords)
  recordlist <- recordlist[n.alignments >= 2]
  unlist(GRangesList(recordlist))
}

.splitread_intersection <- function(x){
  w <- width(intersect(x[1], x[2]))
  ## return 0 if no overlap
  if(length(w) == 0) w <- 0L 
  w
}

splitreadIntersection <- function(g){
  sapply(g, .splitread_intersection)
}

#' Identify rearranged reads -- initiallly unmapped reads that can be
#' aligned by blat to span a novel sequence junction.
#'
#' @return a `GRangesList` of blat records that map to both sides of a
#' sequence  junction. Each list element corresponds to one read that is
#'  aligned to two  locations (i.e., each element of the list consists of the
#'  vector of reads that supports one rearrangement).
#'
#' @export
#'
#' @param linked_bins a \code{GRanges} of linked bins (e.g., gotten by \code{linkedBins(rearrangement.list)})
#'
#' @param blat a data.frame of blat alignment records
#'
#' @param maxgap this maximum gap between the mapped read and the
#'   genomic intervals of the improper read clusters
rearrangedReads <- function(linked_bins, blat, maxgap=500){
  ##lb <- linkedBins(rlist)
  ## filter blat alignments
  is_na <- is.na(blat$Tstart)
  if(any(is_na)){
    blat <- blat[!is_na, ]
  }
  blat <- blat[blat$Tname %in% seqlevels(linked_bins), ]
  blat_gr <- blat_to_granges(blat, linked_bins)
  genome(blat_gr) <- genome(linked_bins)
  ##
  ## A blat record must overlap one of the intervals in a linked bin
  ##
  is.overlap <- overlapsLinkedBins(blat_gr, linked_bins, maxgap=maxgap)
  blat_gr <- blat_gr[ is.overlap ]
  ## We are looking for split read alignments--a read must have at
  ## least 2 alignments.  Get rid of all reads with only a single
  ## alignment, and all alignments that have a match score greater
  ## than 95.
  blat_gr <- multipleAlignmentRecords2(blat_gr)
  ##
  ## Both intervals in a linked bin must overlap a blat record
  ##
  overlaps_both <- overlapsBlatRecord(linked_bins, blat_gr, maxgap)
  sr_list <- vector("list", length(linked_bins))
  names(sr_list) <- names(linked_bins)
  for(i in seq_along(linked_bins)){
    if(!overlaps_both[i]) sr_list[[i]] <- empty_record()
    sr_list[[i]] <-   rearrangedReads2(linked_bins[i], blat_gr, maxgap)
  }
  sr.grl <- GRangesList(sr_list)
  sr.grl
}

empty_record <- function(){
  g <- GRanges()
  g$revmap <- IntegerList()
  g$qStarts <- integer()
  g$blockSizes <- integer()
  g$match <- numeric()
  g$rear.id <- character()
  g$qname <- character()
  g$gapbases <- integer()
  g$Qsize <- integer()
  g
}

blat_to_granges <- function(blat, lb){
  GRanges(blat$Tname, IRanges(blat$Tstart, blat$Tend),
          match=blat$match, qname=blat$Qname,
          qstart=blat$Qstart,
          qend=blat$Qend,
          tStarts=blat$tStarts,
          blockSizes=blat$blockSizes,
          gapbases=blat$Tgapbases,
          blockcount=blat$blockcount,
          qStarts=blat$qStarts,
          Qsize=blat$Qsize,
          seqinfo=seqinfo(lb),
          strand=blat$strand)
}

overlapsLinkedBins <- function(blat_gr, lb, maxgap=500){
  overlapsAny(blat_gr, lb, maxgap=maxgap) |
    overlapsAny(blat_gr, lb$linked.to, maxgap=maxgap)  
}

get_rearrangement_id <- function(records, lb, maxgap=500){
  records$rear.id <- NA
  h1 <- findOverlaps(records, lb, maxgap=maxgap)
  records$rear.id[queryHits(h1)] <- names(lb)[subjectHits(h1)]
  h2 <- findOverlaps(records, lb$linked.to, maxgap=maxgap)
  records$rear.id[queryHits(h2)] <- names(lb)[subjectHits(h2)]
  records$rear.id
}


rearrangedReads2 <- function(linked_bins, blat_gr, maxgap=500){
  lb <- linked_bins
  if(is.null(names(lb))) stop("Expect linked bins to be named by rearrangement id")
  record_overlaps <- overlapsAny(blat_gr, lb, maxgap=maxgap) |
    overlapsAny(blat_gr, lb$linked.to, maxgap=maxgap)
  records <- blat_gr[ record_overlaps ]
  if(length(records) == 0){
    return(GRanges())
  }
  ##
  ## Assign each sequence (qname) to a rearrangement
  ##
  records$rear.id <- get_rearrangement_id(records, lb)
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
  blat.grl <- eachBlockAsGRanges(records_qname)

  # Pass along seqinfo
  seqlevels(blat.grl) <- seqlevels(lb)
  seqlengths(blat.grl) <- seqlengths(lb)
  genome(blat.grl) <- genome(lb)

  overlapFun <- function(g1, lb, ...){
    o1 <- any(overlapsAny(g1, lb, ...))
    o2 <- any(overlapsAny(g1, lb$linked.to, ...))
    o1 && o2
  }
  is.overlap <- sapply(blat.grl, overlapFun, lb=lb, maxgap=500)
  blat.grl <- blat.grl[ is.overlap ]
  blat.grl <- blat.grl [ elementNROWS(blat.grl) == 2 ]
  ##
  ## The split should involve nearly non-overlapping subsequences of
  ## the read, and the total match score should be high (near 100)
  ##
  splitread_ranges <- lapply(blat.grl, sequenceRanges2)
  splitread_int <- splitreadIntersection(splitread_ranges)
  splitread_match <- as.integer(sapply(splitread_ranges, function(x){
    min(sum(x$match), x$Qsize)
  }))
  p <- splitread_match/blat_gr$Qsize[1] * 100
  is_splitread <- splitread_int < 10 & p > 90
  records <- unlist(blat.grl)
  names(records) <- NULL
  ##records_rear <- split(records, records$rear.id)
  ##GRangesList(records_rear)
  records
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
#' \dontrun{
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
#' }
removeAmbiguousAln <- function(rlist, blat_list){
  mapply(removeAmbiguous, x=rlist, blat=blat_list)
}

breakpointInvA <- function(blat.gr, query.gr){
  message("Rearrangement_A")
  ## end of query range is first element
  break1 <- GRanges(seqnames(blat.gr)[1],
                    IRanges(start(blat.gr)[1], width=1),
                    strand=strand(blat.gr[1]))
  break1$linkedto <- GRanges(seqnames(blat.gr)[2],
                             IRanges(start(blat.gr)[2], width=1),
                             strand=strand(blat.gr[2]))
  break1
}

breakpointInvB <- function(blat.gr, query.gr){
  message("Rearrangement_B")
  ## reverse so that lowest query start point is first
  blat.gr <- blat.gr[order(start(query.gr))]
  if(cstrand(blat.gr)[1]=="+"){
    brk1 <- end(blat.gr)[1]
  } else brk1 <- start(blat.gr)[1]
  if(cstrand(blat.gr)[2]=="-"){
    brk2 <- end(blat.gr)[2]
  } else brk2 <- start(blat.gr)[2]
  break1 <- GRanges(seqnames(blat.gr)[1],
                    IRanges(brk1, width=1),
                    strand=strand(blat.gr[1]))
  break1$linkedto <- GRanges(seqnames(blat.gr)[2],
                             IRanges(brk2, width=1),
                             strand=strand(blat.gr[2]))
  break1
}

breakpointsInversion <- function(x, rear){
  lb <- linkedBins(rear)
  blat.gr <- GRanges(x$Tname, IRanges(x$Tstart, x$Tend),
                     strand=x$strand,
                     match=x$match,
                     blockcount=x$blockcount,
                     blockSizes=x$blockSizes,
                     seqinfo=seqinfo(lb),
                     qname=x$Qname)
  blat.gr <- blat.gr[!duplicated(blat.gr)]
  query.gr <- GRanges(x$Tname, IRanges(x$Qstart, x$Qend),
                      strand=x$strand,
                      match=x$match,
                      size=x$Qsize,
                      qname=x$Qname)
  query.gr <- query.gr[!duplicated(query.gr)]

  blat.grl <- split(blat.gr, blat.gr$qname)
  query.grl <- split(query.gr, query.gr$qname)
  brks <- NULL
  for(i in seq_along(blat.grl)){
    blat.gr2 <- blat.grl[[i]]
    query.gr2 <- query.grl[[i]]
    ##
    ## Rearrangement_A proximal part of inversion
    ## - genomic start for second half of query overlaps linkedBins interval
    ##   -- breakpoint is in linkedbins
    ##  (Figure S9, Ordulu)
    ## - genomic start of the first half of query overlaps linkedTo interval
    ##  (Figure S9, Ordulu)
    ix <- order(-end(query.gr2))
    blat.gr2 <- blat.gr2[ix]
    query.gr2 <- query.gr2[ix]
    isRearrangementA <- overlapsAny(blat.gr2[1], lb, maxgap=200)
    if(isRearrangementA) {
      brks2 <- breakpointInvA(blat.gr2, query.gr2)
      brks2$rearrangement <- "A"
      brks <- c(brks, brks2)
    } else{
      ## must be rearrangement_B
      brks2 <- breakpointInvB(blat.gr2, query.gr2)
      brks2$rearrangement <- "B"
      brks <- c(brks, brks2)
    }
    brks
  }
  brks2 <- unlist(GRangesList(brks))
}

breakpointsForQuery <- function(x, seqinfo){
  ## For each split read supporting rearrangement,
  ## determine the sequence junction
  blat.gr <- GRanges(x$Tname, IRanges(x$Tstart, x$Tend),
                     strand=x$strand,
                     match=x$match,
                     blockcount=x$blockcount,
                     blockSizes=x$blockSizes,
                     seqinfo=seqinfo)
  blat.gr <- blat.gr[!duplicated(blat.gr)]
  query.gr <- GRanges(x$Tname, IRanges(x$Qstart, x$Qend),
                      strand=x$strand,
                      match=x$match,
                      size=x$Qsize)
  query.gr <- query.gr[!duplicated(query.gr)]
  ##
  ## order by chromosome and then by alignment with highest score
  ##
  ix <- order(factor(as.character(seqnames(blat.gr)),
                     levels=seqlevels(blat.gr)),
              -blat.gr$match)
  blat.gr2 <- blat.gr[ix]
  if(sum(width(blat.gr2)) > query.gr$size[1]){
    overhang <- sum(width(blat.gr2)) - query.gr$size[1]
    if(cstrand(blat.gr2)[1] == "+"){
      end(blat.gr2)[1] <- end(blat.gr2)[1] - overhang
    } else {
      start(blat.gr2)[1] <- start(blat.gr2)[1] + overhang
    }
  }
  query.gr2 <- query.gr[ix]
  breakpoints <- ifelse(cstrand(blat.gr2)=="+",
                        end(blat.gr2),
                        start(blat.gr2))
  breaks.gr <- GRanges(seqnames(blat.gr2), IRanges(breakpoints, width=1),
                       strand=strand(blat.gr2))
  breaks.gr
}

iscnName <- function(A){
  extdata <- system.file("extdata", package="SNPchip")
  cytobands <- read_tsv(file.path(extdata, "cytoBand_hg19.txt"),
                        col_names=FALSE) %>%
    set_colnames(c("seqnames", "start", "end", "band", "stain")) %$%{
      GRanges(seqnames, IRanges(start, end), band=band, stain=stain)
    }
  band1 <- subsetByOverlaps(cytobands, A) %>%
    .$band
  band2 <- subsetByOverlaps(cytobands, A$linkedto) %>%
    .$band
  ## for A, end of query range is listed in the GRanges
  ##  the linkedto field contains the first part of the query
  if(cstrand(A$linkedto) == "+"){
    if(substr(band2, 1, 1) == "p"){
      terminus1 <- "pter->"
    } else {
      terminus1 <- "pter_cen->"
    }
  } else {
    if(substr(band2, 1, 1) == "p"){
      terminus1 <- "qter_cen->"
    } else {
      terminus1 <- "qter_->"
    }
  }
  if(cstrand(A) == "+"){
    if(substr(band1, 1, 1) == "p"){
      terminus2 <- "->pter"
    } else {
      terminus2 <- "->cen_pter"
    }
  } else {
    if(substr(band1, 1, 1) == "p"){
      terminus2 <- "->_cen_qter"
    } else {
      terminus2 <- "->qter"
    }
  }
  chrA <- gsub("chr", "", chromosome(A))
  iscn <- paste0("seq[hg19]inv(", chrA, ")(",
                 band1, ";", band2, ")")
  ## go from start of read to end of read
  ## for A, end of read is in linkedBins and start of read is in 'linkedto'
  hgsv <- paste0(terminus1, band1,
                 "(", start(A$linkedto), ")::",
                 band2, band1, "(",
                 start(A), "-", start(A$linkedto), ")::",
                 band2, "(", start(A), ")",
                 terminus2)
  A$iscn <- iscn
  A$hgsv <- hgsv
  A
}
