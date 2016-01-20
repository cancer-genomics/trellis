#' @include Deletion-class.R
NULL

#' Accessor for 'copy' assays
#'
#' Extract matrix of log2-transformed estimates of copy number
#' relative to autosomal mode or median
#'
#' @export
#' @docType methods
#' @param object A \code{RangedSummarizedExperiment}
#' @rdname copynumber-methods
setGeneric("copynumber", function(object) standardGeneric("copynumber"))

#' @rdname copynumber-methods
#' @export
setGeneric("copynumber<-", function(object, value) standardGeneric("copynumber<-"))

#' @rdname copynumber-methods
#' @aliases copynumber,RangedSummarizedExperiment-method
setMethod(copynumber, "RangedSummarizedExperiment",
          function(object){
            assays(object)$copy
          })

#' @param value matrix of copy number estimates
#' @rdname copynumber-methods
#' @aliases copynumber<-,RangedSummarizedExperiment,matrix-method
setReplaceMethod("copynumber", signature(object="RangedSummarizedExperiment", value="matrix"),
          function(object, value){
            assays(object)$copy <- value
            object
          })

#' @rdname copynumber-methods
#' @export
#' @aliases copynumber,StructuralVariant-method
setMethod("copynumber", "StructuralVariant", function(object) object@copynumber)

#' @rdname copynumber-methods
#' @export
#' @aliases copynumber<-,StructuralVariant,ANY-method
setReplaceMethod("copynumber", "StructuralVariant", function(object, value){
  object@copynumber <- value
  object
})
