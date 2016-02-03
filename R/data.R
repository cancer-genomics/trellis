#' BLAT alignment of reads that were unmapped by ELAND
#'
#' This is a subset of reads that were involved in an unmapped-mapped
#' pair.  The mapped reads were near a rearrangement supported by
#' improperly paired mapped-mapped reads ('mapped-mapped' indicates
#' both reads were alignmed by ELAND).  We used BLAT to assess whether
#' a split read alignment for the unmapped read supports the new
#' sequence junction.  
#' 
#' @docType data
#' @name blat_unmapped
#' @aliases blat_unmapped
#' @rdname blat_unmapped
#' @format a \code{data.frame} of blat records
NULL

#' An example RearrangementList object
#'
#' A \code{RearrangementList} object contains genomic intervals that
#' are linked by improper read pairs where both reads were mapped.  In
#' addition, it also contains the improper reads.
#'
#' @seealso \code{\link[svclasses]{RearrangementList}}
#' 
#' @docType data
#' @name rearrangement_list
#' @aliases rearrangement_list
#' @rdname rearrangement_list
#' @format a \code{RearrangementList} object
NULL
