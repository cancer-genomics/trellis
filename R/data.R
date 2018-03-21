#' A deletions object of class StructuralVariant
#'
#' A \code{StructuralVariant} object containing data supporting
#' deletions.
#'
#' @docType data
#' @keywords datasets
#' @name deletions
#' @usage data(deletions)
#' @aliases deletions
#'
#' @format a \code{StructuralVariant} object containing several
#'   deletions.
#'
#' @examples
#' data(deletions)
#' variant(deletions)
#' improper(deletions)
#' deletions[1]
#' improper(deletions[1])
#' calls(deletions)
NULL

#' An example AmpliconGraph
#'
#' An \code{AmpliconGraph} for sample CGOV2T.
#'
#' @docType data
#' @keywords datasets
#' @name amplicons-data
#' @usage data(amplicons)
#' @format a \code{AmpliconGraph} object.
#' @examples
#' data(amplicons)
#' ampliconRanges(amplicons)
NULL

#' An example RearrangementList
#'
#' A \code{RearrangementList} object.
#'
#' @docType data
#' @keywords datasets
#' @name rear_list-data
#' @usage data(rear_list)
#' @format a \code{RearrangementList}
#' @aliases rear_list
#' @examples
#' data(rear_list)
#' head(modalRearrangement(rear_list))
#' linkedBins(rear_list)
NULL

#' An example Transcripts object
#'
#' @docType data
#' @keywords datasets
#' @name rear_cds-data
#' @usage data(rear_cds)
#' @aliases rear_cds
#' @format a \code{Transcripts} object
#' @examples
#' data(rear_cds)
NULL

#' An AmpliconGraph object
#'
#' @name amplicon_graph
#' @docType data
NULL

#' A StructuralVariant object
#'
#' @name deletion
#' @docType data
NULL

#' A StructuralVariant object
#'
#' @name deletions2
#' @docType data
NULL

#' A GRanges object of segmentation data
#'
#' @name segments
#' @docType data
NULL

#' A list of preprocess data for calling deletion/amplification
#'
#' @name pdata
#' @docType data
NULL

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
#' @seealso \code{\link{RearrangementList}}
#'
#' @docType data
#' @name rearrangement_list
#' @aliases rearrangement_list
#' @rdname rearrangement_list
#' @format a \code{RearrangementList} object
NULL
