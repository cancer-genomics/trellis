#' @include help.R
NULL

#' Extract sequences of reads supporting rearrangements from a bam file
#'
#'
#' @export
#' @return a data.frame
#' @param object a \code{Rearrangement} or \code{RearrangementList} object
#' @param bam.file complete path to BAM file
#' @param params a \code{RearrangementParams} object
#'
#' @param MAX the maximum number of read pairs to extract for a
#'   rearrangement.  If the number of read pairs supporting a
#'   rearrangement is greater than MAX, a random sample of MAX
#'   supporting read pairs is returned.
#'
#' @rdname getSequenceOfReads-methods
setGeneric("getSequenceOfReads", function(object, bam.file,
                                          params=RearrangementParams(),
                                          MAX=25L, build)
  standardGeneric("getSequenceOfReads"))


# #' Parse BAM file for improper read pairs near a set of GRanges
# #'
# #' All reads aligned to the intervals given by
# #' \code{queryRanges(object)} are identified by the low-level function
# #' \code{.scan_all_readpairs}.  This function reads alignments by
# #' \code{readGAlignments} and then makes pairs of the alignments by
# #' \code{makeGAlignmentPairs2}.  The latter function is an adaption of
# #' the function \code{makeGAlignmentPairs} implemented in the
# #' \code{GenomeAlignments} package but allows for the read pairs to be
# #' improper.
# #'
# #' @param object Typically an \code{AmpliconGraph}, though the only
# #'   requirement is that the method \code{queryRanges} is defined
# #' @param bam.file character-vector providing valid complete path to a
# #'   bam file
# #' @param flags length-two integer vector as given by \code{scanBamFlags}
# #' @export
# setGeneric(name = "get_readpairs", def = function(object, bam.file, flags=scanBamFlag()) {standardGeneric("get_readpairs")})
# setMethod(
#   f = "get_readpairs",
#   signature = "AmpliconGraph",
#   definition = function(object, bam.file, flags=scanBamFlag()) {
#     g <- queryRanges(object)
#     .get_readpairs2(g, bam.file, flags)
#   }
# )
# setMethod(
#   f = "get_readpairs",
#   signature = "GRanges",
#   definition = function(object, bam.file, flags=scanBamFlag()) {
#     .get_readpairs2(g, bam.file, flags)
#   }
# )






