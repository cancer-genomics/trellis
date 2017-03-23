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
                                          MAX=25L)
  standardGeneric("getSequenceOfReads"))
