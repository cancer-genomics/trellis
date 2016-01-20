#' @include help.R
NULL

#' GRangesList-derived class for deletions
#'
#' Each element in the list are the set of deletions for a sample
#' 
#' @export
#' @rdname DeletionList-methods
#' @param object see \code{showMethods("DeletionList")}
setGeneric("DeletionList", function(object) standardGeneric("DeletionList"))

setClass("DeletionList", contains="GRangesList")

#' @rdname DeletionList-methods
#' @aliases DeletionList,list-method
setMethod("DeletionList", "list", function(object){
  DeletionList(GRangesList(object))
})

#' @rdname DeletionList-methods
#' @aliases DeletionList,GRangesList-method
setMethod("DeletionList", "GRangesList", function(object){
  as(object, "DeletionList")
})

#' @rdname DeletionList-methods
#' @aliases DeletionList,missing-method
setMethod("DeletionList", "missing", function(object){
  as(GRangesList(), "DeletionList")
})
