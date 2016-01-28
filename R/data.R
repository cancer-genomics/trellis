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

#' An amplicons object of class AmpliconGraph
#'
#' An \code{AmpliconGraph} of amplicons in sample CGOV2T.
#'
#' @seealso \code{\link[svplots]{plot_amplicons}}
#' @docType data
#' @keywords datasets
#' @name amplicons-data
#' @usage data(amplicons)
#' @format a \code{AmpliconGraph} object.
#' @examples
#' data(amplicons)
#' ampliconRanges(amplicons)
NULL
