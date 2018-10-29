blat_colnames <- c("match",
                   "mis-match",
                   "rep.match",
                   "N's",
                   "Qgapcount",
                   "Qgapbases",
                   "Tgapcount",
                   "Tgapbases",
                   "strand",
                   "Qname",
                   "Qsize",
                   "Qstart",
                   "Qend",
                   "Tname",
                   "Tsize",
                   "Tstart",
                   "Tend",
                   "blockcount",
                   "blockSizes",
                   "qStarts",
                   "tStarts")

#' Read output from the command-line blat
#'
#' The command-line version of blat can be downloaded from sourceforge:
#' \url{http://sourceforge.net/projects/blat/files/}
#'
#' @export
#' @param filename character string providing full path to blat output
#' @param ... additional arguments to \code{read_tsv}
#' @return a \code{data.frame} of alignment records from blat
#'
#' @seealso See \code{\link{blatScores}} for evaluating whether the
#'   blat alignments support a novel sequence junction
#'
readBlat <- function(filename, skip=5, col_names=FALSE, ...){
##  blat <- read.delim(filename, skip=2, nrows=5, stringsAsFactors=FALSE, sep="\t", header=FALSE)
##  nms <- paste0(as.character(blat[1, ]), as.character(blat[2, ]))
##  nms <- gsub(" ", "", nms)
##  ##nms[22:23] <- c("seq1", "seq2")
  blat <- read.delim(filename, skip=5, stringsAsFactors=FALSE,
                     header=FALSE, sep="\t")
  colnames(blat) <- blat_colnames
  ##  colnames(blat) <- nms
  ##  blat <- read_table(filename, skip=skip, col_names=col_names,
  ##                     col_types=col_types, ...) %>%
  ##    set_colnames(blat_colnames) %>%
  ##    mutate(Tname=gsub(".fa", "", Tname))
  ##blat2 <- as.tibble(blat)
  blat
}

blatGRanges <- function(blat, sl=paste0("chr", c(1:22, "X", "Y", "M"))){
  g <- GRanges(blat$Tname, IRanges(blat$Tstart, blat$Tend),
               strand=factor(blat$strand, levels=c("+", "-", "*")))
  seqlevels(g, pruning.mode="coarse") <- sl
  g
}

blatGRanges2 <- function(blat, sl=paste0("chr", c(1:22, "X", "Y", "M"))){
  chr <- blat$Tname
  starts <- blat$Tstart
  ends <- blat$Tend
  strand <- factor(blat$strand, levels=c("+", "-", "*"))
  g <- GRanges(chr, IRanges(starts, ends),
               blockcount=blat$blockcount,
               Qname=blat$Qname,
               rid=blat$rid,
               score=blat$score)
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

#' Annotate blat records
#'
#' @param blat \code{data.frame} of blat records from \code{readBlat}
#' @param tag.sequences \code{data.frame} of read sequences
#'
#' @examples
#' data(tags)
#' extdata <- system.file("extdata", package="svbams")
#' blat.file <- file.path(extdata, "blat_alignment.txt")
#' blat_aln <- readBlat(blat.file)
#' records <- annotateBlatRecords(blat_aln, tags)
#' annotateBlatRecords(records, tags)
#' 
#' @export
annotateBlatRecords <- function(blat, tag.sequences){
  Qname <- NULL
  qname2 <- paste0(tag.sequences$qname, "_", tag.sequences$read)
  query.sequences <- setNames(tag.sequences$seq, qname2)
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
  ## Tsize from blat is the length of the chromosome.  Change to be the length of the reference sequence that matches the query sequence (read)
  blat$Tsize <- blat$Tend-blat$Tstart
  blat$is_overlap <- genomeAndBlatOverlap(blat)

  Tstart <- Tend <- Qsize <- rearrangement.id <- NULL
  tags <- tag.sequences %>%
    unite("Qname", c("qname", "read"))
  rid <- tags %>%
    group_by(Qname) %>%
    summarize(rid=unique(rearrangement.id))
  blat2 <- left_join(blat, rid, by="Qname") %>%
    mutate(score=match/Qsize)
  rownames(blat2) <- NULL
  blat2 <- as.tibble(blat2)
  blat2
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

##blatStatsPerTag <- function(blat.records, tag_length){
##  stats <- blat.records %>%
##    mutate(is_90 = match > 90,
##           Tsize=abs(Tend-Tstart),
##           is_size_near100=Tsize > (tag_length-1/5*tag_length) &
##             Tsize < (tag_length + 1/5*tag_length),
##           is_90 = is_90 & is_size_near100) %>%
##    group_by(Qname) %>%
##    summarize(n.matches=sum(is_90),
##              n.eland.matches=sum(is_90 & is_overlap)) %>%
##    mutate(is_pass=n.matches == 1 & n.eland.matches==1)
##  return(stats)
##  if(FALSE){
##    ## summary statistics for a single read
##    is_90 <- blat.records$match > 90
##    ## Tsize in blat.records is size of target sequence (chromosome)
##    Tsize <- abs(blat.records$Tend-blat.records$Tstart)
##    is_size_near100 <- Tsize > (tag_length - 1/5*tag_length) &
##      Tsize < (tag_length + 1/5*tag_length)
##    is_90 <- is_90 & is_size_near100
##    is_overlap <- blat.records$is_overlap
##    ## calculate number of near-perfect matches for each tag
##    n.matches <- sapply(split(is_90, blat.records$Qname), sum)
##    n.eland.matches <- sapply(split(is_overlap & is_90, blat.records$Qname), sum)
##    ## there should only be one eland alignment with high quality (n.matches == 1)
##    ## AND this read should overlap the eland alignment
##    n.matches==1 & n.eland.matches==1
##  }
##}

##.blatStatsRearrangement <- function(blat, thr=0.8, tag_length){
##  cols <- c("Qname", "match", "is_overlap", "Tstart", "Tend")
##  blat <- blat[, cols]
##  if(nrow(blat) == 0) return(NULL)
##  splitby <- blat$Qname %>%
##    strsplit("R[12]") %>%
##    sapply("[", 1)
##  blat$tag_index <- blat$Qname %>% strsplit(splitby) %>%
##    sapply("[", 2)
##  ##blat$tag_index <- addXCoordinateForTag(blat)
##  ##stats <- blatStatsPerTag(blat, tag_length)
##  Tend <- Tstart <- Tsize <- is_90 <- is_size_near100 <- NULL
##  is_overlap <- n.matches <- n.eland.matches <- NULL
##  stats <- blat %>%
##    mutate(is_90 = match > 90,
##           Tsize=abs(Tend-Tstart),
##           is_size_near100=Tsize > (tag_length-1/5*tag_length) &
##             Tsize < (tag_length + 1/5*tag_length),
##           is_90 = is_90 & is_size_near100) %>%
##    group_by(Qname) %>%
##    summarize(n.matches=sum(is_90),
##              n.eland.matches=sum(is_90 & is_overlap)) %>%
##    mutate(is_pass=n.matches == 1 & n.eland.matches==1)
##
##  proportion_pass <- mean(stats)
##  qnames <- gsub("R[12]", "", blat$Qname)
##  n_tags <- length(unique(qnames))
##  is_pass <- proportion_pass > thr & n_tags >= 5
##  ##blat$passQC <- rep(proportion_pass > thr, nrow(blat))
##  blat$passQC <- is_pass
##  blat
##}



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
#' data(tags)
#' extdata <- system.file("extdata", package="svbams")
#' blat.file <- file.path(extdata, "blat_alignment.txt")
#' blat_aln <- readBlat(blat.file)
#' blat <- annotateBlatRecords(blat_aln, tags)
#' blatScores(blat, tags)
#' @export
#' @param blat a \code{data.frame} of results from  command-line blat
#' @param tags a \code{data.frame} containing read names and the original alignment locations
#' @param prop.pass a length-one numeric vector indicating the fraction of
#'   reads at a rearrangement that must pass the read-level QC.
#' @param min.tags the minimum number of tags that pass BLAT QC for each rearrangement
#' @param id sample id
blatScores <- function(blat, tags, id, min.tags=5, prop.pass=0.8){
  blockcount <- tsize <- Qsize <- Qname <- overlaps_genome <- NULL
  rid <- number_alignments <- p_overlap_genome <- number_tags <- NULL
  Tend <- Tstart <- Tsize <- tag_length <- is_size_near100 <- NULL
  is_90 <- n.matches <- n.genome.matches <- is_pass <- proportion_pass <- ngats <- NULL
  tags <- tags %>%
    unite("Qname", c("qname", "read"))
  filter <- dplyr::filter
  blat2 <- blat %>%
    filter(score >= 0.90 & blockcount==1) %>%
    filter(Tsize >= (Qsize - 1/5*Qsize) & Tsize <= (Qsize + 1/5*Qsize))
  blat.g <- blatGRanges2(blat2)
  ##
  ## Record for each blat record whether the blat alignment corresponds to the whole-genome-aligner
  ##
  genome.g <- GRanges(tags$seqnames, IRanges(tags$start, tags$end),
                      Qname=tags$Qname, rid=tags$rearrangement.id)
  ##
  ## Force to have same seqlevels
  ##
  seqlevels(blat.g) <- seqlevels(genome.g) <- paste0("chr", c(1:22, "X", "Y", "M"))
  hits <- findOverlaps(blat.g, genome.g, ignore.strand=TRUE, maxgap=500)
  keep <- blat.g$Qname[queryHits(hits)] == genome.g$Qname[subjectHits(hits)]
  hits <- hits[keep]
  blat2$overlaps_genome <- FALSE
  blat2$overlaps_genome[ queryHits(hits) ] <- TRUE
  ##
  ## We've already filtered low scoring alignments
  ## - number below is number of high scoring alignments
  ##
  blat3 <- blat2 %>%
    group_by(Qname)  %>%
    summarize(number_alignments=n(),
              overlaps_genome=sum(overlaps_genome),
              rid=unique(rid)) %>%
    ungroup %>%
    group_by(rid) %>%
    summarize(##number_reads=length(unique(Qname)),
      number_tags=n(),
      p_overlap_genome=mean(overlaps_genome/number_alignments)) %>%
    mutate(passQC=p_overlap_genome >= prop.pass &
             number_tags >= min.tags) %>%
    mutate(rid=factor(rid))
  blat3
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
  if(is.numeric(x)) return(x)
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
  ## We are looking for split read alignments--a read with a blockcount of 1
  ## must have at least 2 alignments. Get rid of all reads with only a single
  ## alignment having a block count of 1. Among these reads, remove alignments that have a match score greater than 95.
  ##
  ## Alternatively, blat can report a single alignment with multiple blocks.
  ## Here, the high match score should be high, two of the blocks should hit the linked bins, and there should be more than start site in the target genome.
  ##
  tstarts <- integer_vector(blat_gr$tStarts)
  is_candidate <- (blat.gr$match < 95 & blat.gr$blockcount == 1) |
    (blat.gr$match > 95 & blat.gr$blockcount > 1 & length(tstarts) > 1)
  ##blat.gr$match < 95 | blat.gr$blockcount > 1
  ##browser()
  ##blat.gr$match < 95 & blat.gr$blockcount == 1
  is_candidate
}

numberAlignmentRecords <- function(blat.gr){
  L <- length(blat.gr)
  if(L == 1){
    L <- blat.gr$blockcount
  }
  L
}

.each_block_granges <- function(g){
  starts <- integer_vector(g$tStarts)
  L <- length(starts)
  widths <- integer_vector(g$blockSizes)
  qstarts <- integer_vector(g$qStarts)
  bsizes <- integer_vector(g$blockSizes)
  chrom <- rep(chromosome(g), g$blockcount)
  qends <- qstarts+bsizes
  qends <- tmp2$Qsize
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
  recordlist <- recordlist[ n.alignments >= 2 ]
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

#' Compute the intersection of the split read alignment
#'
#' Example
#'
#' split 1     :  -----------
#' split 2     :           -----------
#' intersection:  2
intersectionAligned <- function(g){
  grl <- split(g, g$blockcount)
  w <- rep(NA, length(grl))
  for(i in seq_along(grl)){
    g2 <- grl[[i]]
    if(g2$blockcount[1] == 1){
      g2 <- IRanges(g2$qstart, g2$qend)
      tmp <- width(intersect(g2[1], g2[2]))
      ## no intersection should be recorded as zero
      w[i] <- ifelse(length(tmp) == 0, 0, tmp)
    } else {
      ## if the read is aligned in blocks and the blocks are non-overlapping, the intersection is zero
      tstarts <- integer_vector(g2$tStarts)
      tends <- tstarts + integer_vector(g2$blockSizes)
      ir <- IRanges(tstarts, tends)
      hits <- findOverlaps(ir, ir)
      hits <- hits[queryHits(hits) != subjectHits(hits)]
      ir1 <- ir[queryHits(hits)]
      ir2 <- ir[subjectHits(hits)]
      widths <- rep(NA, length(hits))
      for(j in seq_along(ir1)){
        widths[j] <- width(intersect(ir1[j], ir2[j]))
      }
      w[i] <- max(widths)
    }
  }
  w
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
  ## BLAT fields
  ## matches - Number of matching bases that aren't repeats.
  ## misMatches - Number of bases that don't match.
  ## repMatches - Number of matching bases that are part of repeats.
  ## nCount - Number of 'N' bases.
  ## qNumInsert - Number of inserts in query.
  ## qBaseInsert - Number of bases inserted into query.
  ## tNumInsert - Number of inserts in target.
  ## tBaseInsert - Number of bases inserted into target.
  ## strand - defined as + (forward) or - (reverse) for query strand. In mouse, a second '+' or '-' indecates genomic strand.
  ## qName - Query sequence name.
  ## qSize - Query sequence size.
  ## qStart - Alignment start position in query.
  ## qEnd - Alignment end position in query.
  ## tName - Target sequence name.
  ## tSize - Target sequence size.
  ## tStart - Alignment start position in target.
  ## tEnd - Alignment end position in target.
  ## blockCount - Number of blocks in the alignment.
  ## blockSizes - Comma-separated list of sizes of each block.
  ## qStarts - Comma-separated list of start position of each block in query.
  ## tStarts - Comma-separated list of start position of each block in target.
  ## filter blat alignments
  is_na <- is.na(blat$Tstart)
  if(any(is_na)){
    blat <- blat[!is_na, ]
  }
  blat <- blat[blat$Tname %in% seqlevels(linked_bins), ]
  blat_gr <- blat_to_granges(blat, seqinfo(linked_bins))
  genome(blat_gr) <- genome(linked_bins)
  ##
  ## A blat record must overlap one of the intervals in a linked bin
  ##
  is.overlap <- overlapsLinkedBins(blat_gr, linked_bins, maxgap=maxgap)
  blat_gr <- blat_gr[ is.overlap ]
  if(length(blat_gr) == 0){
    message("No blat records overlap the linked regions")
    return(NULL)
  }
  ## We are looking for split read alignments--a read must have at
  ## least 2 alignments.  Get rid of all reads with only a single
  ## alignment, and all alignments that have a match score greater
  ## than 95.
  blat_gr <- multipleAlignmentRecords2(blat_gr)
  ##
  ## Total of splits should be near Qsize
  ##  - intersection of splits should be no more than 10% of Qsize
  ##
  grl <- split(blat_gr, names(blat_gr))
  MIN.SIZE <- ceiling(0.95*blat_gr$Qsize[1])
  MAX.INTERSECT <- floor(0.1*blat_gr$Qsize[1])
  queryAligned <- function(g) sum(g$qend - g$qstart)
  total_size_aligned <- sapply(grl, queryAligned)
  intersection_split <- sapply(grl, intersectionAligned)
  size.intersection.filter <- total_size_aligned >= MIN.SIZE &
    intersection_split <= MAX.INTERSECT
  blat_gr <- unlist(grl[ size.intersection.filter ])
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

blat_to_granges <- function(blat, seqinfo){
  GRanges(blat$Tname, IRanges(blat$Tstart, blat$Tend),
          match=blat$match,
          qname=blat$Qname,
          qstart=blat$Qstart,
          qend=blat$Qend,
          tStarts=blat$tStarts,
          blockSizes=blat$blockSizes,
          gapbases=blat$Tgapbases,
          blockcount=blat$blockcount,
          qStarts=blat$qStarts,
          Qsize=blat$Qsize,
          seqinfo=seqinfo,
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
  ## Each query aligned with a single block should have a unique rearrangement id
  ##
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
    ## a query only map to lb or the linked range
    ## lb and the linked range must have at least one hit
    o1 <- overlapsAny(g1, lb, ...)
    o2 <- overlapsAny(g1, lb$linked.to, ...)
    is.overlap <- any((o1 & !o2) | (o2 & !o1)) &&
      sum(o1) >= 1 && sum(o2) >= 1
    is.overlap
  }
  is.overlap <- sapply(blat.grl, overlapFun, lb=lb, maxgap=500)
  blat.grl <- blat.grl[ is.overlap ]
  if(length(blat.grl) == 0){
    msg <- "Split reads do not uniquely align to both ends of the target sequence"
    message(msg)
    return(NULL)
  }
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
  . <- band <- stain <- NULL
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
