#' @include Deletion-class.R
NULL

setGeneric("first<-", function(x, value) standardGeneric("first<-"))

setGeneric("last<-", function(x, value) standardGeneric("last<-"))


##--------------------------------------------------
##
## constructors
##
##--------------------------------------------------

#' Constructor for an empty GAlignmentPairs object
#' 
#' @return an empty \code{GAlignmentPairs} object
#' @export
#' @param first a \code{GAlignments}
#' @param last a \code{GAlignments}
#' @param isProperPair logical for whether to require that the reads
#'   be in a proper pair
#' @rdname GAlignmentPairs
.GAlignmentPairs <- function() {
  GAlignmentPairs(first=GAlignments(),
                  last=GAlignments(),
                  isProperPair=logical())
}


#' Represent read pairs as segments
#'
#' When both reads of a pair are aligned to the same chromosome, this
#' function makes a \code{GRanges} object of the interval (DNA
#' fragment).
#'
#' 
#' @return a \code{GRanges} object
#' @export
#' @param object a \code{GAlignmentPairs} object
readPairsAsSegments <- function(object){
  intrachrom <- chromosome(first(object)) == chromosome(last(object))
  firstIsLeft <- start(first(object)) < start(last(object))
  starts <- start(first(object))
  starts[!firstIsLeft] <- start(last(object))[!firstIsLeft]
  ends <- end(last(object))
  ends[!firstIsLeft] <- end(first(object))[!firstIsLeft]
  all(starts < ends)
  rpsegs <- GRanges(chromosome(first(object))[intrachrom],
                    IRanges(starts[intrachrom],
                            ends[intrachrom]))
  rpsegs
}

setReplaceMethod("first", "GAlignmentPairs", function(x, value){
  x@first <- value
  x
})

setReplaceMethod("last", "GAlignmentPairs", function(x, value){
  x@last <- value
  x
})
