#' A class for storing segmentation parameters
#'
#' @slot alpha A length-one numeric vector
#' @slot undo.splits A length-one character vector
#' @slot undo.SD A length-one numeric vector
#' @slot verbose A length-one numeric vector
#' @export
setClass("SegmentParam", representation(alpha="numeric",
                                        undo.splits="character",
                                        undo.SD="numeric",
                                        verbose="numeric"))

#' A parameters class for segmentation by DNAcopy
#'
#' @seealso \code{\link[DNAcopy]{segment}}
#' 
#' A description for each of these argments can be found by reading
#'   the documentation for the \code{segment} function in the
#'   \code{DNAcopy} package.
#'
#' @examples
#' sp <- SegmentParam()
#' cbs_alpha(sp)
#' undo.splits(sp)
#' undo.SD(sp)
#'
#' @return \code{SegmentParam} object
#' @param alpha a length-one numeric vector 
#' @param undo.splits character string
#' @param undo.SD length-one numeric vector
#' @param verbose length-one numeric vector
#' @rdname SegmentParam-class
#' @export
SegmentParam <- function(alpha=0.001,
                         undo.splits="sdundo",
                         undo.SD=2,
                         verbose=0) {
  new("SegmentParam", alpha=alpha, undo.splits=undo.splits,
      undo.SD=undo.SD, verbose=verbose)
}

#' @param x \code{SegmentParam} object
#' @rdname SegmentParam-class
#' @export
setGeneric("cbs_alpha", function(x) standardGeneric("cbs_alpha"))

#' @rdname SegmentParam-class
#' @export 
setGeneric("undo.splits", function(x) standardGeneric("undo.splits"))

#' @rdname SegmentParam-class
#' @export 
setGeneric("undo.SD", function(x) standardGeneric("undo.SD"))

#' @aliases cbs_alpha,SegmentParam-method
#' @rdname SegmentParam-class
setMethod(cbs_alpha, "SegmentParam", function(x) x@alpha)

#' @aliases undo.splits,SegmentParam-method
#' @rdname SegmentParam-class
setMethod(undo.splits, "SegmentParam", function(x) x@undo.splits)

#' @aliases undo.SD,SegmentParam-method
#' @rdname SegmentParam-class
setMethod(undo.SD, "SegmentParam", function(x) x@undo.SD)
