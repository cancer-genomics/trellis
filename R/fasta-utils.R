.writefasta <- function(x, fasta.file){
  snms <- x$snms
  ##cat("", file=fasta.file)
  for(k in seq_along(snms)){
    cat(">", snms[k], "\n", append=TRUE, sep="", file=fasta.file)
    cat(x$seq[[k]], "\n", append=TRUE, sep="", file=fasta.file)
  }
}

#' @export
writeUnmappedToFasta <- function(tags, fasta.file){
  unlink(fasta.file)
  .writefasta(tags, fasta.file)
}


#' @export
writeToFasta <- function(tags, fasta.file){
  unlink(fasta.file)
  tag.list <- seqsByRearrangement(tags)
  lapply(tag.list, .writefasta, fasta.file=fasta.file)
  file.exists(fasta.file)
}

seqsByRearrangement <- function(x){
  if(is.null(x)) return(x)
  x$snms <- paste(x$qname, x$read, sep="_")
  x <- x[!duplicated(x$snms), ]
  xlist <- split(x, factor(x$rearrangement.id, levels=unique(x$rearrangement.id)))
  xlist <- lapply(xlist, function(x) x[order(x$qname), ])
  xlist
}


fastaFiles <- function(dirs, ids){
  path <- file.path(dirs[["fasta"]], "fasta")
  file.path(path, paste0(ids, ".fa"))
}


# Extract unmapped reads with a mapped mate from a bam file
#
# This function extracts all unmapped reads with a mate that overlaps
# with a set of query genomic intervals.  Internally, this function
# uses \code{scanBam} to scan the bam files.  In our application, we
# typically create a bam file containing only mapped-unmapped read
# pairs as this greatly reduces the size of the bam file to query.
# In particular, we create a bam file with the following set of flags:
#
#       \code{samtools view -b -f 4 -F 8 $input > "unmapped-mapped/${input}"}
#
# @details
#
# The \code{GRanges} object returned by this function includes the
#   sequence of the reads so that the sequences can be subsequently
#   written to disk in fasta format and realigned with a local
#   alignment algorithm such as \code{BLAT} that allows for split
#   read alignments.
#
# @return a \code{GRanges} object of mapped reads with unmapped
#   mates.
#
# @param bam.file character-string providing complete path to BAM file
#
# @param query a \code{GRanges} representation of genomic intervals
#   to query for a mapped read with an unmapped mate.  For example, a
#   set of rearrangement intervals.
#
# @param yield_size the number of reads to extract from the bam file at once using \code{scanBam}
#
# @param maxgap the gap allowed between the query interval and the
#   mapped read to consider the two intervals overlapping
#
# @param what a character vector of fields to keep from the bam file.
#   Defaults to \code{scanBamWhat()}.
#
# @seealso See \code{\link[Rsamtools]{scanBam}} for scanning a bam
#   file for reads matching a set of flags.  See
#   \code{fasta_unmapped} for writting these sequences to disk in
#   fasta format.
# unmapped_read <- function(bam.file, query, yield_size=1e6, maxgap=200,
#                          what=scanBamWhat()){
#   ##
#   ## if more than one bam file, just use the first
#   ##
#   flags <- scanBamFlag(isUnmappedQuery=TRUE, hasUnmappedMate=FALSE, isDuplicate=FALSE) 
#   param <- ScanBamParam(flag=flags, what=what)
#   bamfile <- BamFile(bam.file, yieldSize=yield_size)
#   i <- 1
#   mate_grl <- list()
#   open(bamfile)
#   while(length((chunk0 <- scanBam(bamfile, param=param)[[1]])$qname)){
#     cat(".")
#     mate_gr <- GRanges(chunk0$mrnm, IRanges(chunk0$mpos, width=100))
#     if(any(overlapsAny(mate_gr, query, maxgap=maxgap))){
#     ##,-----------------------------------------------------------
#     ##| Extract sequences to blat. Because 'query' is the unmapped
#     ##| read, I think 'seq' is the sequence of the query (?)
#     ##| 
#     ##`-----------------------------------------------------------
#     mate_gr$seq <- chunk0$seq
#     names(mate_gr) <- chunk0$qname
#     mate_gr <- subsetByOverlaps(mate_gr, query, maxgap=maxgap)
#     mate_grl[[i]] <- mate_gr
#     i <- i+1
#     }
#   }
#   close(bamfile)
#   mate_gr <- unlist(GRangesList(mate_grl))
#   if(FALSE){
#     si <- seqinfo(bamRanges(bview))
#     mate_gr <- keepSeqlevels(mate_gr, seqlevels(si))
#     seqinfo(mate_gr) <- si
#   }
#   mate_gr$snms <- names(mate_gr)
#   mate_gr$seq <- as.character(mate_gr$seq)
#   mate_gr
# }


#' Extract unmapped reads with a mapped mate from a bam file
#'
#' This function extracts all unmapped reads with a mate that overlaps
#' with a set of query genomic intervals.  Internally, this function
#' uses \code{scanBam} to scan the bam files.  In our application, we
#' typically create a bam file containing only mapped-unmapped read
#' pairs as this greatly reduces the size of the bam file to query.
#' In particular, we create a bam file with the following set of flags:
#'
#'       \code{samtools view -b -f 4 -F 8 $input > "unmapped-mapped/${input}"}
#'
#' @details
#'
#' The \code{GRanges} object returned by this function includes the
#'   sequence of the reads so that the sequences can be subsequently
#'   written to disk in fasta format and realigned with a local
#'   alignment algorithm such as \code{BLAT} that allows for split
#'   read alignments.
#'
#' @return a \code{GRanges} object of mapped reads with unmapped
#'   mates.
#'
#' @param bam.file character-string providing complete path to BAM file
#'
#' @param query a \code{GRanges} representation of genomic intervals
#'   to query for a mapped read with an unmapped mate.  For example, a
#'   set of rearrangement intervals.
#'
#' @param yield_size the number of reads to extract from the bam file at once using \code{scanBam}
#'
#' @param maxgap the gap allowed between the query interval and the
#'   mapped read to consider the two intervals overlapping
#'
#' @param what a character vector of fields to keep from the bam file.
#'   Defaults to \code{scanBamWhat()}.
#'
#' @seealso See \code{\link[Rsamtools]{scanBam}} for scanning a bam
#'   file for reads matching a set of flags.  See
#'   \code{fasta_unmapped} for writting these sequences to disk in
#'   fasta format.
#' @examples 
#' extdata <- system.file("extdata", package="svbams")
#' bam <- file.path(extdata, "cgov44t_revised.bam")
#' region <- GRanges(seqnames = "chr8",
#' ranges = IRanges(start = 128691748, end = 128692097))
#' unmapped_read(bam.file = bam, query = region, yield_size = 1e6)
#' @export
unmapped_read <- function(bam.file, query, yield_size=1e6, maxgap=200, what=scanBamWhat()){
  flags <- scanBamFlag(isUnmappedQuery=TRUE, hasUnmappedMate=FALSE, isDuplicate=FALSE) 
  param <- ScanBamParam(flag=flags, what=what, which = query+500)
  bamfile <- BamFile(bam.file, yieldSize=yield_size)
  mate_grl <- list()
  open(bamfile)
  while(length(unlist(sapply(chunk0 <- scanBam(bamfile, param=param), function(x) x[[1]])))) {
    cat(".")
    for (i in 1:length(chunk0)) {
      mate_gr <- GRanges(chunk0[[i]]$mrnm, IRanges(chunk0[[i]]$mpos, width=100))
      if(any(overlapsAny(mate_gr, query, maxgap=maxgap))){ # This should be the case -- can remove if statement
        mate_gr$seq <- chunk0[[i]]$seq
        names(mate_gr) <- chunk0[[i]]$qname
        mate_gr <- subsetByOverlaps(mate_gr, query, maxgap=maxgap)
        mate_grl[[length(mate_grl)+1]] <- mate_gr
      } 
    }
  }
  close(bamfile)
  mate_gr <- unlist(GRangesList(mate_grl))
  if(FALSE){
    si <- seqinfo(bamRanges(bview))
    mate_gr <- keepSeqlevels(mate_gr, seqlevels(si))
    seqinfo(mate_gr) <- si
  }
  mate_gr$snms <- names(mate_gr)
  mate_gr$seq <- as.character(mate_gr$seq)
  mate_gr
}
















#' Write the read sequences of unmapped-mapped read pairs to disk in fasta format
#'
#' @seealso See \code{\link{unmapped_read}} for extracting a
#'   \code{GRanges} representation of unmapped-mapped read pairs where
#'   the mapped read overlap with a set of query genomic intervals,
#'   such as candiate rearrangements.
#'
#' @return nothing is returned
#' 
#' @export
#' @param dp a \code{DataPaths} object
#' @param unmapped_read a \code{GRanges} object
#' @param id a length-one character vector for the sample id
fasta_unmapped <- function(dp, unmapped_read, id){
  message("Writing unmapped-mapped read sequences to fasta files")
  outfile <- file.path(dp[["fasta_unmapped"]], paste0(id, ".fa"))
  if(file.exists(outfile)) return(invisible())
  .writefasta(unmapped_read, outfile)
  return(invisible())
}
