#' @include Rearrangement-class.R
NULL

#' @aliases improper,Rearrangement-method
#' @rdname improper-methods
setMethod("improper", "Rearrangement", function(object) object@improper)

#' @rdname StructuralVariant-class
#' @aliases improper,StructuralVariant-method
setMethod("improper", "StructuralVariant", function(object){
  i <- unlist(indexImproper(object))
  object@improper[i]
})

#' @rdname improper-methods
#' @aliases improper,StructuralVariant,GAlignmentPairs-method
#' @param value a \code{GAlignmentPairs} object
#' @keywords internal
setReplaceMethod("improper", c("StructuralVariant","GAlignmentPairs"),
                 function(object,value){
                   object@improper <- value
                   index <- initializeImproperIndex2(variant(object), value)
                   indexImproper(object) <- index
                   object
                 })

#' @rdname improper-methods
#' @aliases improper,Rearrangement,GAlignmentPairs-method
#' @param value a \code{GAlignmentPairs} object
setReplaceMethod("improper", c("Rearrangement", "GAlignmentPairs"),
                 function(object, value) {
                   object@improper <- value
                   object
                 })
