#' @include Rearrangement-class.R
NULL

#' Accessor for improper read pairs
#'
#' Extract the improperly paired reads.  Reads are flagged as improper
#' by the alignment algorithm.  Typically, this indicates that the
#' separation between a read and its mate is larger than expected, the
#' orientation of the read and its mate is different from that
#' anticipated in the reference genome, or the reads align to
#' different chromosomes.
#' 
#' 
#' @export
#' @rdname improper-methods
#' @seealso \code{\linkS4class{Rearrangement}}
#' @param object a \code{Rearrangement} or \code{StructuralVariant} object
setGeneric("improper", function(object) standardGeneric("improper"))

#' @rdname improper-methods
#' @export
setGeneric("improper<-", function(object, value) standardGeneric("improper<-"))

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
