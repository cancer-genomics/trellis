#' @include AllGenerics.R
NULL

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
