#' @include Deletion-class.R
NULL

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
