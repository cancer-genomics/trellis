#' @include Rearrangement-class.R
NULL

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

#' A container for a list of Rearrangement objects
#'
#' A given sample may have multiple rearrangements.  It is convenient
#' to encapsulate the data supporting each rearrangement in a
#' \code{Rearrangement} object, and to list all of these
#' rearrangements in a \code{RearrangementList}.
#'
#' @examples
#' RearrangementList()
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

#' @rdname fractionLinkingTags
#' @aliases fractionLinkingTags,RearrangementList-method
setMethod("fractionLinkingTags", "RearrangementList", function(object){
  sapply(object, fractionLinkingTags)
})

#' @aliases length,RearrangementList-method
#' @rdname RearrangementList-class
setMethod("length", "RearrangementList", function(x) length(x@data))

#' @aliases linkedBins,RearrangementList-method
#' @rdname linkedBins-methods
setMethod("linkedBins", "RearrangementList", function(object){
  dat <- object@data
  lblist <- lapply(dat, linkedBins)
  nms <- sapply(lblist, names)
  grl <- GRangesList(lblist)
  g <- unlist(grl)
  names(g) <- nms
  g
})

#' @aliases linkedBins,RearrangementList,ANY-method
#' @rdname linkedBins-methods
setReplaceMethod("linkedBins", "RearrangementList", function(object, value){
  dat <- object@data
  for(i in seq_along(dat)){
    linkedBins(dat[[i]]) <- value[i]
  }
  object@data <- dat
  object
})

#' @rdname modalRearrangement
#' @aliases modalRearrangement,RearrangementList-method
setMethod("modalRearrangement", "RearrangementList", function(object){
  sapply(object, modalRearrangement)
})


#' @aliases names,RearrangementList-method
#' @rdname RearrangementList-class
setMethod("names", "RearrangementList", function(x) x@names)

#' @aliases numberLinkingRP,RearrangementList-method
#' @rdname numberLinkingRP-methods
#' @export
setMethod("numberLinkingRP", "RearrangementList", function(object){
  sapply(object, numberLinkingRP)
})

#' @rdname RearrangementList-class
#' @aliases overlapsAny,RearrangementList,GRanges-method
#' 
#' @param query a \code{RearrangementList}
#' @param subject a \code{GRanges} object
#' @param ... additional arguments ignored
setMethod("overlapsAny", c("RearrangementList", "GRanges"), function(query, subject, ...){
  lb <- linkedBins(query)
  overlaps_firstbin <- overlapsAny(lb, subject)
  overlaps_secondbin <- overlapsAny(lb$linked.to, subject)
  overlaps_firstbin | overlaps_secondbin
})

#' @rdname percentRearrangement
#' @aliases percentRearrangement,RearrangementList-method
setMethod("percentRearrangement", "RearrangementList", function(object){
  sapply(object, percentRearrangement)
})

#' @aliases sapply,RearrangementList-method
#' @rdname sapply-methods
setMethod("sapply", "RearrangementList", function(X, FUN, ..., simplify=TRUE, USE.NAMES=TRUE){
  results <- setNames(rep(NA, length(X)), names(X))
  for(i in seq_along(X)){
    results[i] <- FUN(X[[i]], ...)
  }
  results
})

#' @aliases lapply,RearrangementList-method
#' @rdname sapply-methods
setMethod("lapply", "RearrangementList", function(X, FUN, ...){
  ##results <- setNames(rep(NA, length(X)), names(X))
  results <- vector("list", length(X))
  setNames(results, names(X))
  for(i in seq_along(X)){
    results[[i]] <- FUN(X[[i]], ...)
  }
  results
})

#' @rdname splitReads
#' @aliases splitReads,RearrangementList,GRangesList-method
setReplaceMethod("splitReads", c("RearrangementList", "GRangesList"),
                 function(object, value){
                   orig_order <- names(object)
                   object2 <- object[ names(object) %in% names(value) ]
                   object2 <- object2 [ names(value) ]
                   for (i in 1:length(object2)) {
                     splitReads(object2[[i]]) <- value[[i]]
                   }
                   notchanged <- object [ !names(object) %in% names(object2) ]
                   object3 <- c(notchanged, object2)
                   object3 <- object3[ orig_order ]
                   object3
                 })


#' @rdname splitReads
#' @aliases splitReads,RearrangementList-method
setMethod("splitReads", "RearrangementList", 
          function(object){
            split_reads <- vector("list", length(object))
            for(i in seq_along(object)){
              split_reads[[i]] <- splitReads(object[[i]])
            }
            grl <- GRangesList(split_reads)
            names(grl) <- names(object)
            grl
          })

##--------------------------------------------------
##
## Subsetting
##
##--------------------------------------------------

#' @rdname RearrangementList-class
#' @aliases $,RearrangementList,ANY-method
setReplaceMethod("$", "RearrangementList", function(x, name, value){
  x@colData[[name]] <- value
  x
})

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

##--------------------------------------------------
##
## combining
##
##--------------------------------------------------
setMethod("c", "RearrangementList",
          function(x, ..., recursive=FALSE){
            args <- list(x, ...)
            if(all(lengths(args) == 0)) return(RearrangementList())
            rear.list <- lapply(args, function(x) x@data)
            rear.list2 <- do.call("c", rear.list)
            nms <- names(rear.list2)
            coldat.list <- lapply(args, colData)
            ## this does not work for some reason
            ##coldat <- do.call(rbind, coldat.list)
            coldat.list <- lapply(coldat.list, as.data.frame)
            coldat <- as(do.call(rbind, coldat.list), "DataFrame")
            new("RearrangementList",
                data=rear.list2,
                names=nms,
                elementType="Rearrangement",
                colData=coldat)
          })
