#' @include Rearrangement-class.R
NULL

#' The StructuralVariant class stores data pertaining to somatic deletions
#'
#' \code{StructuralVariant} is a vector-like class for storing genomic
#' intervals of somatic deletions and relevant supporting data.
#' Properly and improperly paired reads and indices so that subsetting
#' methods (\code{[}) extracts the relevant data for a given deletion.
#'
#' @slot variant A  \code{GRanges} object
#' @slot proper A  \code{GAlignmentPairs} object of properly paired reads
#' @slot improper A \code{GAlignmentPairs} object of improperly paired
#'   reads
#'
#' @slot copynumber A \code{numeric} vector of the mean log ratios
#'   corresponding to the \code{variant} intervals. The numeric vector
#'   must have the same length as \code{variant}.
#'
#' @slot calls A length-N \code{character} vector of deletion calls
#'   where N is the number of intervals in \code{variant}.  Possible
#'   calls are (i) hemizygous: hemizygous deletion with fewer than 5
#'   supporting improper read pairs (ii) hemizygous+: hemizygous
#'   deletion with 5 or more supporting improper read pairs (iii)
#'   homozygous: homozygous deletion with fewer than 5 supporting
#'   improper read pairs (iv) homozygous+: homozygous deletion with 5
#'   or more supporting improper read pairs
#'
#' @slot index_proper a list of indices.  Each element is a vector of
#'   integers.  The kth element are the elements of the
#'   \code{GAlignmentPairs} object in slot \code{proper} supporting
#'   the kth deletion.
#'
#' @slot index_improper A list of indices.  Each element is a vector
#'   of integers.  The kth element are the elements of the
#'   \code{GAlignmentPairs} object in slot \code{improper} supporting
#'   the kth deletion.
#'
#' @slot grouped_variant a factor indicating whether deletions belong
#'   to the same group (e.g., overlapping hemizygous deletions)
#'
#' @slot length_improper Keeps track of the numer of the total number
#'   of improper read pairs for quick summaries of this object.
#'
#' @slot length_proper Keeps track of the numer of the total number of
#'   proper read pairs for quick summaries of this object.
#'
#' @examples
#' StructuralVariant()
#'
#'
#' @export
setClass("StructuralVariant",
         representation(variant="GRanges",
                        proper="GAlignmentPairs",
                        improper="GAlignmentPairs",
                        copynumber="numeric",
                        calls="character",
                        index_proper="list",
                        index_improper="list",
                        grouped_variant="factor",
                        length_improper="integer",
                        length_proper="integer"))

max_proper_read_index <- function(object){
  i <- sapply(indexProper(object), function(x) {
    is_null <- is.null(x)
    if(is_null) return(0)
    max(x)
  })
  i
}

is_valid_proper_index <- function(object){
  msg <- TRUE
  if(sum(elementNROWS(indexProper(object))) == 0){
    return(msg)
  }
  i <- max_proper_read_index(object)
  maxi <- max(i)
  if(maxi > length(object@proper)){
    msg <- "Out of bounds indexing for proper read pairs"
  }
  msg
}

setValidity("StructuralVariant", function(object){
  msg <- TRUE
  if(length(variant(object))==0){
    return(msg)
  }
  L <- length(object)
  L1 <- length(indexProper(object))
  L2 <- length(indexImproper(object))
  L3 <- length(copynumber(object))
  L4 <- length(calls(object))
  if(L1 != L || L2 != L) {
    msg <- "The index list for proper and improper pairs must be the same as the number of SVs"
    return(msg)
  }
  if(L3 != L || L4 != L){
    msg <- paste("The copynumber vector and calls vector must be the same length\n",
                 "as the number of variants")
    return(msg)
  }
  id <- names(variant(object))
  id2 <- names(indexProper(object))
  id3 <- names(indexImproper(object))
  if(!identical(id, id2) || !identical(id, id3)){
    msg <- paste("Names of index list for proper and improper pairs\n",
                 "must be the same as the names for the SVs")
    return(msg)
  }
  if(sum(elementNROWS(indexImproper(object))) > 0){
    ##i <- as.integer(unlist(indexImproper(object)))
    i <- sapply(indexImproper(object), function(x) {
      is_null <- is.null(x)
      if(is_null) return(0)
      max(x)
    })
    maxi <- max(i)
    if(maxi > object@length_improper){
      msg <- "Out of bounds indexing for improper read pairs"
      return(msg)
    }
  }
  msg <- is_valid_proper_index(object)
  msg
})

#' @aliases [,StructuralVariant,numeric,ANY-method
#' @rdname StructuralVariant-class
setMethod("[", signature(x="StructuralVariant", i="numeric"),
          function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    x@variant <- variant(x)[i]
    x@calls <- calls(x)[i]
    x@copynumber <- copynumber(x)[i]
    x@index_improper <- indexImproper(x)[i]
    x@index_proper <- indexProper(x)[i]
    x@grouped_variant <- x@grouped_variant[i]
  }
  x
})

setMethod("show", "StructuralVariant", function(object){
  cat("StructuralVariant class\n")
  cat("    #SVs         :", length(object), "\n")
  cat("    proper RPs   :", sum(elementNROWS(indexProper(object))), "\n")
  cat("    improper RPs :", sum(elementNROWS(indexImproper(object))), "\n")
  ##cat("    proper RPs   :", object@length_proper, "\n")
  ##cat("    improper RPs :", object@length_improper, "\n")
})

##
## Constructor
##

#' @export
#' @param variant a \code{GRanges} object of the deletion genomic intervals
#' @param proper a \code{GAlignmentPairs} object of proper read pairs near the deletion
#' @param improper a \code{GAlignmentPairs} object of improper read pairs near the deletion
#' @param copynumber the mean log2 ratio (derived from read depth) of the deletion region
#' 
#' @param calls character vector indicating the type of deletion
#'   (hemizygous, hemizygous+, homozygous , homozygous+, overlapping
#'   hemizygous)
#' 
#' @param index_proper integer vector used internally for subsetting
#' @param index_improper integer vector used internally for subsetting
#' @param grouped_variant a grouping factor for the deletion intervals
#' @rdname StructuralVariant-class
StructuralVariant <- function(variant=GRanges(),
                              proper=.GAlignmentPairs(),
                              improper=.GAlignmentPairs(),
                              copynumber,
                              calls,
                              index_proper,
                              index_improper,
                              grouped_variant){
  nms <- names(variant)
  if(is.null(nms) && length(variant) > 0){
    nms <- paste0("sv", seq_len(length(variant)))
    names(variant) <- nms
  } ##else nms <- character()
  L <- length(variant)
  if(missing(index_proper))
    index_proper <- setNames(vector("list", L), nms)
  if(missing(index_improper))
    index_improper <- setNames(vector("list", L), nms)
  if(missing(grouped_variant))
    grouped_variant <- factor(rep(NA, L))
  if(missing(copynumber)) copynumber <- as.numeric(rep(NA, L))
  if(missing(calls)) calls <- as.character(rep(NA, L))
  length_improper <- length(improper)
  length_proper <- length(proper)
  new("StructuralVariant",
      variant=variant,
      proper=proper,
      improper=improper,
      copynumber=copynumber,
      calls=calls,
      index_proper=index_proper,
      index_improper=index_improper,
      grouped_variant=grouped_variant,
      length_improper=length_improper,
      length_proper=length_proper)
}

##--------------------------------------------------
##
## methods
##
##--------------------------------------------------


#' @rdname calls-method
#' @aliases calls,StructuralVariant-method
#' @export
setMethod("calls", "StructuralVariant", function(object) object@calls)

#' @rdname calls-method
#' @aliases calls,StructuralVariant,ANY-method
#' @export
setReplaceMethod("calls", "StructuralVariant", function(object, value){
  object@calls <- value
  object
})

initializeImproperIndex2 <- function(gr, improper_rp, param=DeletionParam()){
  hits <- findOverlaps(gr, improper_rp, minimumGapWidth(param))
  subj_hits <- subjectHits(hits)
  index_improper <- split(subj_hits, names(gr)[queryHits(hits)])
  result <- setNames(vector("list", length(gr)), names(gr))
  result[names(index_improper)] <- index_improper
  result
}

.match_index_variant <- function(index_improper, gr, value){
  i <- match(names(value), names(gr))
  index_improper[i] <- value
  index_improper
}

#' @aliases groupedVariant,StructuralVariant-method
#' @rdname groupedVariant-method
setMethod("groupedVariant", "StructuralVariant", function(object) object@grouped_variant)

#' @aliases groupedVariant,StructuralVariant,ANY-method
#' @rdname groupedVariant-method
setReplaceMethod("groupedVariant", "StructuralVariant", function(object, value) {
  object@grouped_variant <- value
  object
})

#' @aliases indexProper,StructuralVariant,list-method
#' @rdname indexing-methods
#' @export
setReplaceMethod("indexProper", c("StructuralVariant","list"),
                 function(object, value){
                   object@index_proper <- value
                   object
##
##                   ## only update the SV indices for which proper pairs are available
##                   index <- match(names(value), names(variant(object)))
##                   object@index_proper[index] <- value
##                   ## indices without any matches assign null
##                   index < which(! names(value) %in% names(variant(object)))
##                   object@index <- NULL
##                   object
                 })

#' @aliases indexImproper,StructuralVariant,list-method
#' @rdname indexing-methods
#'@export
setReplaceMethod("indexImproper", c("StructuralVariant","list"),
                 function(object, value){
                   ##                    index_improper <- .match_index_variant(indexImproper(object),
                   ##                                                           variant(object),
                   ##                                                           value)
                   ##                    object@index_improper <- index_improper
                   object@index_improper <- value
                   object
                 })


#' @rdname indexing-methods
#' @aliases indexProper,StructuralVariant-method
#' @keywords internal
setMethod("indexProper", "StructuralVariant", function(object) object@index_proper)


#' @rdname indexing-methods
#' @aliases indexProper,StructuralVariant,list-method
#' @keywords internal
setReplaceMethod("indexProper", c("StructuralVariant","list"),
                 function(object, value){
                   ## only update the SV indices for which proper pairs are available
                   index <- match(names(value), names(variant(object)))
                   object@index_proper[index] <- value
                   object
                 })

#' @rdname indexing-methods
#' @aliases indexImproper,StructuralVariant-method
#' @keywords internal
setMethod("indexImproper", "StructuralVariant", function(object) object@index_improper)

.match_index_variant <- function(index_improper, gr, value){
  i <- match(names(value), names(gr))
  index_improper[i] <- value
  index_improper
}

#' @rdname indexing-methods
#' @aliases indexImproper,StructuralVariant,list-method
#' @keywords internal
setReplaceMethod("indexImproper", c("StructuralVariant","list"),
                 function(object, value){
                   index_improper <- .match_index_variant(indexImproper(object),
                                                          variant(object),
                                                          value)
                   object@index_improper <- index_improper
                   object
                 })


#' @rdname StructuralVariant-class
#' @aliases length,StructuralVariant-method
setMethod("length", "StructuralVariant", function(x) length(variant(x)))

#' @rdname StructuralVariant-class
#' @aliases names,StructuralVariant-method
#' @param x a \code{StructuralVariant}
setMethod("names", "StructuralVariant", function(x) names(variant(x)))

#' Extract number of improper reads for each element in a StructuralVariant
#'
#'
#' @param object a \code{StructuralVariant}
#' @export
#' @examples
#' data(deletions)
#' numberImproper(deletions)
numberImproper <- function(object){
  num <- rep(NA, length(object))
  lengths(indexImproper(object))
}

#' @rdname proper-methods
#' @aliases proper,StructuralVariant-method
setMethod("proper", "StructuralVariant", function(object){
  i <- unlist(indexProper(object))
  object@proper[i]
})

#' @rdname proper-methods
#' @aliases proper,StructuralVariant,GAlignmentPairs-method
setReplaceMethod("proper", c("StructuralVariant","GAlignmentPairs"),
                 function(object,value){
                   object@proper <- value
                   object
                 })

#' @rdname StructuralVariant-class
#' @aliases readPairs,StructuralVariant-method
setMethod("readPairs", "StructuralVariant", function(object){
  prp <- proper(object)
  irp <- improper(object)
  if(!identical(colnames(mcols(first(prp))),
                colnames(mcols(first(irp))))){
    mcols(first(prp)) <- NULL
    mcols(first(irp)) <- NULL
    mcols(last(prp)) <- NULL
    mcols(last(irp)) <- NULL    
  }
  rp <- c(prp, irp)
  i <- order(start(first(rp)))
  rp <- rp[i]
  rp
})



#' @rdname StructuralVariant-class
#' @export
#' @aliases sort,StructuralVariant-method
#' @param decreasing logical
#' @param ... ignored
setMethod("sort", "StructuralVariant", function(x, decreasing=FALSE, ...){
  v <- variant(x)
  i <- order(v)
  x <- StructuralVariant(variant=v[i],
                         calls=calls(x)[i],
                         copynumber=copynumber(x)[i],
                         proper=x@proper,
                         improper=x@improper,
                         index_proper=indexProper(x)[i],
                         index_improper=indexImproper(x)[i],
                         grouped_variant=groupedVariant(x)[i])
  x
})


#' @rdname StructuralVariant-class
#' @export
#' @aliases variant,StructuralVariant-method
#' @param object a \code{StructuralVariant} 
setMethod("variant", "StructuralVariant", function(object) object@variant)

#' @rdname StructuralVariant-class
#' @export
#' @aliases variant<-,StructuralVariant,ANY-method
#' @keywords internal
setReplaceMethod("variant", "StructuralVariant", function(object, value){
  object@variant <- value
  object
})

##--------------------------------------------------
##
## Subsetting
##
##--------------------------------------------------

#' @rdname StructuralVariant-class
#' @aliases [,StructuralVariant,ANY,ANY-method
setMethod("[", "StructuralVariant", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    x@variant <- variant(x)[i]
    x@calls <- calls(x)[i]
    x@copynumber <- copynumber(x)[i]
    x@index_improper <- indexImproper(x)[i]
    x@index_proper <- indexProper(x)[i]
    x@grouped_variant <- x@grouped_variant[i]
  }
  x
})

#' @rdname StructuralVariant-class
#' @aliases [[,StructuralVariant,ANY,ANY-method
setMethod("[[", "StructuralVariant", function(x, i){
  if(!missing(i)){
    x@variant <- variant(x)[i]
    x@calls <- calls(x)[i]
    x@copynumber <- copynumber(x)[i]
    x@index_improper <- indexImproper(x)[i]
    x@index_proper <- indexProper(x)[i]
    x@grouped_variant <- x@grouped_variant[i]
  }
  x
})

##
## Vector like operations
##

#' @rdname sapply-methods
#' @aliases lapply,StructuralVariant-method
#' @export
setMethod("lapply", "StructuralVariant", function(X, FUN, ...){
  FUN <- match.fun(FUN)
  result <- vector("list", length(X))
  for(i in seq_along(X)){
    result[[i]] <- FUN(X[i], ...)
  }
  result
})

#' Functionals for structural variant classes
#'
#' @keywords internal 
#' @rdname sapply-methods
#' @aliases sapply,StructuralVariant-method
#' @param X see \code{showMethods("sapply")}
#' @param FUN a function
#' @param ... additional arguments to \code{FUN}
#' @param simplify logical
#' @param USE.NAMES logical
#' @export
setMethod("sapply", "StructuralVariant", function(X, FUN, ..., simplify=TRUE, USE.NAMES=TRUE){
  FUN <- match.fun(FUN)
  answer <- lapply(X = X, FUN = FUN, ...)
  names(answer) <- names(X)
  if (USE.NAMES && is.character(X) && is.null(names(answer)))
    names(answer) <- X
  if (!identical(simplify, FALSE) && length(answer))
    simplify2array(answer, higher = (simplify == "array"))
  else answer
})

setMethod("rename", "StructuralVariant", function(x, ...){
  nms <- paste0("sv", seq_along(variant(x)))
  names(variant(x)) <- nms
  ip <- indexProper(x)
  names(ip) <- nms
  x@index_proper <- ip
  ip <- indexImproper(x)
  names(ip) <- nms
  x@index_improper <- ip
  names(copynumber(x)) <- nms
  x
})

#' @rdname StructuralVariant-class
#' @aliases variant<-,StructuralVariant,GRanges-method
setReplaceMethod("variant", c("StructuralVariant", "GRanges"),
                 function(object, value){
                   object@variant <- value
                   object
                 })

setMethod("combine", c("StructuralVariant", "StructuralVariant"),
          function(x, y, ...){
            cnv <- c(variant(x), variant(y))
            prp <- c(proper(x), proper(y))
            irp <- c(improper(x), improper(y))
            cn <- c(copynumber(x), copynumber(y))
            calls <- c(calls(x), calls(y))

            prp_index <- initializeProperIndex2(cnv, prp, zoom.out=1)
            irp_index1 <- initializeImproperIndex2(cnv, irp, DeletionParam())
            irp_index2 <- updateImproperIndex2(cnv, irp, maxgap=2e3)
            irp_index3 <- .match_index_variant(irp_index1, cnv, irp_index2)
            sv <- StructuralVariant(variant=cnv,
                                    proper=prp,
                                    improper=irp,
                                    copynumber=cn,
                                    calls=calls,
                                    index_proper=prp_index,
                                    index_improper=irp_index3)
            sv
          })

setMethod("reviseJunction", "StructuralVariant", function(object){
  v <- variant(object)
  if(length(v) > 1) stop("method only defined for length-one StructuralVariant")
  irp <- improper(object)
  if(length(irp) < 5) return(variant(object))
  irp2 <- leftAlwaysFirst(irp)
  left_border <- end(first(irp2))
  right_border <- start(last(irp2))
  d <- right_border-left_border
  if(any(d > 500)){
    st <- max(left_border[d > 500])+50
    en <- min(right_border[d > 500])-50
    if(st > en) {
      st <- max(left_border)
      en <- min(right_border)
      if(st > en){
        ##stop("Interval for SV is too small after revision -- start
        ##location is larger than end location")
        ##
        ## Keep original boundary
        ##
        st <- start(v)
        en <- end(v)
      }
    }
  } else return(v)
  v <- GRanges(seqnames(v), IRanges(st, en))
  v
})
