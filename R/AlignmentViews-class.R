#' @include Preprocess-class.R
NULL

#' Accessor for file names containing improper read pair alignments
#'
#' @return a character-vector of file paths
#' @rdname AlignmentViews2-class
#' @export
#' @param object a \code{AlignmentViews2} object
setGeneric("improperPaths", function(object) standardGeneric("improperPaths"))

setGeneric("improperPaths<-", function(object, value) standardGeneric("improperPaths<-"))
