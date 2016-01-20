#' @include help.R
NULL

#' A container for a list of Rearrangement objects
#'
#' A given sample may have multiple rearrangements.  It is convenient
#' to encapsulate the data supporting each rearrangement in a
#' \code{Rearrangement} object, and to list all of these
#' rearrangements in a \code{RearrangementList}.
#' 
#' @export
#' @slot data a list of \code{Rearrangement} objects
#' @slot elementType a character vector
#' @slot names a charcter vector of names.  Must have the same length as the list object.
#' @slot colData a \code{DataFrame}
#' @slot modal_rearrangement a \code{character} vector
#' @slot percent_rearrangement a \code{numeric} vector
setClass("RearrangementList", representation(data="list",
                                             elementType="character",
                                             names="character",
                                             colData="DataFrame",
                                             modal_rearrangement="character",
                                             percent_rearrangement="numeric"))


setValidity("RearrangementList", function(object){
  msg <- TRUE
  elements <- sapply(object@data, class)
  if(!all(elements==elementType(object))){
    msg <- "All elements of list must be of type Rearrangement"
    return(msg)
  }
  msg
})

setMethod("show", "RearrangementList", function(object){
  cat("An object of class 'RearrangementList'\n")
  cat("  number rearrangement objects: ", length(object), "\n")
  cat("  Use '[[i]]' to return a single Rearrangement object'\n")
})


#' Constructor for \code{RearrangementList} class
#'
#' @return a \code{RearrangementList} object
#' @rdname RearrangementList-class
#' @export
#' @keywords internal
#' @param object a list of \code{Rearrangement} objects
#' 
#' @param modal_rearrangement a character vector of the modal
#'   rearrangement corresponding to each element in \code{object}.
#'   Must be the same length as \code{object}.
#' 
#' @param percent_rearrangement a numeric vector indicating the
#'   fraction of improper read pairs supporting the modal
#'   rearrangement type.  Must be the same length as \code{object}.
#' 
#' @param colData a \code{DataFrame} containing metadata on the
#'   \code{Rearrangement} objects.
#' 
setGeneric("RearrangementList", function(object,
                                         modal_rearrangement,
                                         percent_rearrangement,
                                         colData) standardGeneric("RearrangementList"))

#' @rdname RearrangementList-class
#' @aliases RearrangementList,missing-method
setMethod("RearrangementList", "missing",
          function(object, colData){
            new("RearrangementList",
                data=list(), names=character(),
                elementType="Rearrangement",
                colData=DataFrame())
          })

#' @rdname RearrangementList-class
#' @aliases RearrangementList,list-method
setMethod("RearrangementList", "list",
          function(object, colData){
            if(missing(colData)) colData <- DataFrame(row.names=names(object))
            new("RearrangementList",
                data=object,
                names=names(object),
                elementType="Rearrangement",
                colData=colData)
          })

listRearrangements2 <- function(object){
  xlist <- setNames(vector("list", length(object)), names(linkedBins(object)))
  for(i in seq_along(object)){
    xlist[[i]] <- object[i]
  }
  rlist <- RearrangementList(xlist)
  rlist
}

#' @rdname RearrangementList-class
#' @aliases RearrangementList,Rearrangement-method
setMethod("RearrangementList", "Rearrangement", function(object) {
  listRearrangements2(object)
})

#' @rdname RearrangementList-class
#' @aliases colData,RearrangementList-method
setMethod("colData", "RearrangementList", function(x, ...){
  x@colData
})

#' @rdname RearrangementList-class
#' @aliases colData,RearrangementList,ANY-method
setReplaceMethod("colData", "RearrangementList", function(x, value){
  x@colData <- value
  x
})

#' @rdname RearrangementList-class
#' @export
#' @keywords internal
setMethod("elementType", "RearrangementList", function(x, ...) x@elementType)

#' @aliases length,RearrangementList-method
#' @rdname RearrangementList-class
setMethod("length", "RearrangementList", function(x) length(x@data))

#' @aliases names,RearrangementList-method
#' @rdname RearrangementList-class
setMethod("names", "RearrangementList", function(x) x@names)

##--------------------------------------------------
##
## Subsetting
##
##--------------------------------------------------

#' @rdname RearrangementList-class
#' @aliases $,RearrangementList-method
setMethod("$", "RearrangementList", function(x, name){
  colData(x)[[name]]
})

#' @return \code{RearrangementList}
#' @rdname RearrangementList-class
#' @param x a \code{RearrangementList} object
#' @param i a numeric or character vector of rearrangement names
#' @param j ignored
#' @param ... ignored
#' @param drop ignored
setMethod("[", "RearrangementList", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    if(is(i, "character")){
      i <- match(i, names(x))
    }
    x@data <- x@data[i]
    x@names <- x@names[i]
    x@modal_rearrangement <- x@modal_rearrangement[i]
    x@percent_rearrangement <- x@percent_rearrangement[i]
    colData(x) <- colData(x)[i, , drop=FALSE]
  }
  x
})

#' @seealso \code{\linkS4class{Rearrangement}}
#' @return a \code{Rearrangement} object
#' @rdname RearrangementList-class
setMethod("[[", "RearrangementList", function(x, i, j, ..., drop=FALSE){
  ## returns a Rearrangement object
  if(!missing(i)){
    if(is(i, "character")){
      i <- match(i, names(x))
    }
    x <- x@data[[i]]
  }
  x
})

#' @seealso \code{\linkS4class{Rearrangement}}
#' @return a \code{Rearrangement} object
#' @rdname RearrangementList-class
setReplaceMethod("[[", "RearrangementList", function(x, i, j, ..., drop=FALSE, value){
  ## returns a Rearrangement object
  if(!missing(i)){
    x@data[[i]] <- value
  }
  x
})
