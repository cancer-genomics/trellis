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
#' modalRearrangement(rear_list)
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
