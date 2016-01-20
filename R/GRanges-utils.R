.expand <- function(granges, size){
  if(any(is.na(seqlengths(granges)))) stop("seqlengths must be specified")
  if(identical(length(size), 1L)){
    size <- rep(size, length(granges))
  }
  starts <- pmax(start(granges)-size, 1L)
  ends <- pmin(end(granges)+size, seqlengths(granges)[chromosome(granges)])
  d1 <- start(granges)-starts
  d2 <- ends-end(granges)
  i <- which(d1 > 0 & d2 > 0 & d1 != d2)
  if(length(i) > 0){
    dd <- pmin(d1[i], d2[i])
    starts[i] <- start(granges)[i]-dd
    ends[i] <- end(granges)[i]+dd
  }
  start(granges) <- starts
  end(granges) <- ends
  granges
}

#' Expand a genomic interval
#'
#' Expand a genomic interval by a user-specified number of basepairs.
#'
#' @examples
#' library(GenomeInfoDb)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' library(GenomicRanges)
#' si <- keepSeqlevels(seqinfo(BSgenome.Hsapiens.UCSC.hg19), "chr1")
#' gr <- GRanges("chr1", IRanges(15e3, 20e3), seqinfo=si)
#' gr2 <- expandGRanges(gr, 5e3)
#' ## will not go beyond the seqlengths
#' gr3 <- GRanges("chr1", IRanges(seqlengths(gr)-5000, seqlengths(gr)-2000),
#'                seqinfo=si)
#' gr4 <- expandGRanges(gr3, 10e3)
#' end(gr4) == seqlengths(gr) 
#' 
#' 
#' @export
#' @param granges a \code{GRanges} object
#' 
#' @param size a length-one numeric vector specifying the number of
#'   basepairs to subtract from the start and add to the ends of the
#'   intervals in object \code{granges}
expandGRanges <- function(granges, size) .expand(granges, size)

#' Find the fraction of genomic interval overlapping with another
#' genomic interval.
#'
#' 
#' Computing the fraction of a query genomic interval overlapping with
#' another set of genomic interval amounts to finding the intersection
#' of the query with the other intervals (subject) and then dividing
#' by the width of the query.
#'
#' @examples
#'  library(GenomicRanges)
#'  autosomes <- function() paste0("chr", 1:22)
#'  del <- GRanges(c("chr1", "chr2"), IRanges(c(1, 1), c(10, 10)))
#'  filters <- GRanges(c("chr1", "chr2"), IRanges(c(5, 8), c(20, 10000)))
#'  test <- intOverWidth(del, filters)
#'  identical(test, c(0.6, 0.3))
#'
#' @return A numeric vector of the same length as the query
#'   \code{GRanges}
#' 
#' @export
#' @param query a \code{GRanges} object
#' @param subject a \code{GRanges} object
intOverWidth <- function(query, subject){
  subject <- reduce(subject)
  p <- rep(0, length(query))
  hits <- findOverlaps(query, subject)
  if(length(hits) == 0) return(p)
  ## vectorized intersect
  i <- queryHits(hits)
  j <- subjectHits(hits)
  intersection <- width(pintersect(query[i], subject[j]))
  numer <- sapply(split(intersection, i), sum)
  denom <- width(query)[as.integer(names(numer))]
  p[as.integer(names(numer))] <- numer/denom
  p
}
