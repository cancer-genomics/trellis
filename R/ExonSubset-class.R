#' @include help.R
NULL

#' @export
#' @keywords-internal
#' @rdname ExonSubset-class
#' @param object a \code{ExonSubset} object
#' @return a \code{LogicalList}
setGeneric("inRearrangement.left", function(object) standardGeneric("inRearrangement.left"))

#' @export
#' @keywords-internal
#' @rdname ExonSubset-class
#' @return a \code{LogicalList}
setGeneric("inRearrangement.right", function(object) standardGeneric("inRearrangement.right"))

#' @export
#' @keywords-internal
#' @rdname ExonSubset-class
setGeneric("name.left", function(object) standardGeneric("name.left"))

#' @export
#' @keywords-internal
#' @rdname ExonSubset-class
setGeneric("name.right", function(object) standardGeneric("name.right"))

#' Container for storing information to subset transcripts
#' 
#' @examples
#' ExonSubset()
#' @export
#' @slot name.left character string
#' @slot inRearrangement.left a \code{LogicalList}
#' @slot name.right character string
#' @slot inRearrangement.right a \code{LogicalList}
setClass("ExonSubset", representation(name.left="character", ##e.g, "A"
                                      inRearrangement.left="LogicalList",
                                      name.right="character",
                                      inRearrangement.right="LogicalList"))

#' @rdname ExonSubset-class
#' @export
#' @param name.left name of left transcript
#' @param inRearrangement.left \code{LogicalList} subsetting the left transcripts
#' @param name.right name of right transcript
#' @param inRearrangement.right \code{LogicalList} subsetting the right transcripts
ExonSubset <- function(name.left=character(),
                       inRearrangement.left=LogicalList(),
                       name.right=character(),
                       inRearrangement.right=LogicalList()){
  new("ExonSubset", name.left=name.left,
      inRearrangement.left=inRearrangement.left,
      name.right=name.right,
      inRearrangement.right=inRearrangement.right)
}

#' @aliases name.left,ExonSubset-method
#' @rdname ExonSubset-class
setMethod("name.left", "ExonSubset", function(object) object@name.left)

#' @aliases name.right, ExonSubset-method
#' @rdname ExonSubset-class
setMethod("name.right", "ExonSubset", function(object) object@name.right)

#' @aliases names,ExonSubset-method
#' @param x an \code{ExonSubset} object
#' @rdname ExonSubset-class
setMethod("names", "ExonSubset", function(x) paste0(name.left(x), name.right(x)))

#' @aliases inRearrangement.left,ExonSubset-method
#' @rdname ExonSubset-class
setMethod("inRearrangement.left", "ExonSubset", function(object) object@inRearrangement.left)

#' @aliases inRearrangement.right,ExonSubset-method
#' @rdname ExonSubset-class
setMethod("inRearrangement.right", "ExonSubset", function(object) object@inRearrangement.right)

setMethod("isEmpty", "ExonSubset", function(x){
  length(name.left(x)) == 0
})

setValidity("ExonSubset", function(object){
  ##
  msg1 <- "The logical list for the 5' protein in a chimeric protein should always be TRUE followed by FALSE, or all FALSE."
  msg2 <- "The logical list for the 2nd transcript in a chimeric protein should always be FALSE followed by TRUE, or all TRUE"
  msg <- TRUE
  lt <- inRearrangement.left(object)
  if(length(lt) > 0){
    is.valid <- vector("list", length(lt))
    ##is.valid <- foreach(x=lt, .combine=c) %do% {
    for(i in seq_along(lt)){
      x <- lt[[i]]
      if(any(x) & any(!x)){
        index.true <- max(which(x))
        index.false <- min(which(!x))
        index.true < index.false
        is.valid[[i]] <- index.true
      } else is.valid[[i]] <- TRUE
    }
    is.valid <- unlist(is.valid)
    if(!all(is.valid)){
      return(msg1)
    }
  }
  rt <- inRearrangement.right(object)
  if(length(rt) > 0){
    is.valid <- vector("list", length(rt))
    ##is.valid <- foreach(x=rt, .combine=c) %do% {
    for(i in seq_along(rt)){
      x <- rt[[i]]
      if(any(x) & any(!x)){
        index.true <- min(which(x))
        index.false <- max(which(!x))
        is.true <- index.true > index.false
        is.valid[[i]] <- is.true
      } else is.valid[[i]] <- TRUE
    }
    is.valid <- unlist(is.valid)
    if(!all(is.valid)){
      return(msg2)
    }
  }
  msg
})
