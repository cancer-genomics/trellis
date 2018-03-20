#' @include AllGenerics.R
NULL

#' Set numeric scale for \code{PreprocessViews2} object
#'
#' Exported only for internal use by other packages.
#' 
#' @return numeric
#' 
#' @export
#' @docType methods
#' @param x A \code{PreprocessViews2} object
#' @param value a length-one numeric vector
#' @rdname setScale-method
#' @keywords internal
setGeneric("setScale<-", function(x, value) standardGeneric("setScale<-"))

#' Indexes the GRanges of a BamViews-derived class
#'
#' @docType methods
#' @rdname indexRanges-method
#' @param object a \code{BamViews}-derived object
#' @export
#' @keywords internal
setGeneric("indexRanges", function(object) standardGeneric("indexRanges"))

#' @param value an integer-vector 
#' @rdname indexRanges-method
#' @export
#' @keywords internal
setGeneric("indexRanges<-", function(object, value) standardGeneric("indexRanges<-"))

#' Accessor for file paths
#'
#' File path to intermediate data saved to disk
#' 
#' @export
#' @param object a \code{PreprocessViews2} object
#' @rdname paths
#' @aliases paths paths<-
setGeneric("paths", function(object) standardGeneric("paths"))

#' @param value a character-vector of file paths to intermediate data
#' @rdname paths
#' @aliases paths,PreprocessViews2-method  paths<-,PreprocessViews2-method
#' @export
setGeneric("paths<-", function(object,value) standardGeneric("paths<-"))

#' A container for storing views of preprocessed data
#'
#' This class directly extends the \code{BamViews} class and provides
#' views to intermediate data (saved to disk) from the analysis of one
#' or more bam files.
#'
#' @seealso See \code{\link{paths}} for the character-vector of file
#'   paths to the intermediate data.
#' 
#' @slot scale a length-one numeric vector.  We scale numeric data by
#'   the value of \code{scale}, round to the nearest integer, and then
#'   save as an integer.  This slot is for internal use only.
#'
#' @examples
#' PreprocessViews2()
#' paths(PreprocessViews2())
#' pviews <- PreprocessViews2()
#' paths(pviews) <- character()
#' paths(pviews)
#'
#' @export
setClass("PreprocessViews2", representation(scale="numeric"),
         contains="BamViews")

setValidity("PreprocessViews2", function(object){
  msg <- TRUE
  if(!"range_index" %in% colnames(mcols(rowRanges(object)))){
    msg <- "mcols of the rowRanges is missing a required column 'range_index'"
    return(msg)
  }
  if(!identical(length(indexRanges(object)), length(rowRanges(object)))){
    msg <- "length of indexRanges and rowRanges must be the same"
    return(msg)
  }
  msg
})

#' Constructor for PreprocessViews2
#'
#' @return A \code{PreprocessViews2} object
#'
#' @export
#' @param object can be \code{missing} or an existing \code{BamViews} object
#' @rdname PreprocessViews2-class
PreprocessViews2 <- function(object){ ## A BamViews object
  if(missing(object)){
    object <- BamViews()
  }
  object <- as(object, "PreprocessViews2")
  indexRanges(object) <- seq_len(nrow(object))
  setScale(object) <- 1 ## by default
  object
}

##
##  Accessors / methods for PreprocessViews2
##

#' Accessor for reading assay data saved to disk
#'
#' We have adapted the assays method to read assay data from disk
#' using a \code{PreprocessViews2} object.
#'
#' 
#' REFACTOR: both ... and withDimnames are ignored.  we should just
#' provide a different generic/method.
#'
#' 
#' @return a R x C matrix, where R is the length of
#'   \code{rowRanges(x)} and C is the number of samples given by
#'   \code{ncol(x)}
#' 
#' @export
#' 
#' @param x a \code{PreprocessViews2} object
#' @param ... ignored
#' @param withDimnames ignored
setMethod("assays", "PreprocessViews2", function(x, ..., withDimnames=FALSE){
  result <- matrix(NA, nrow(x), ncol(x))
  for(j in seq_len(ncol(x))){
    result[, j] <- readRDS(paths(x)[j])[indexRanges(x)]
  }
  colnames(result) <- colnames(x)
  if(getScale(x) != 1){
    result <- result/getScale(x)
  }
  result
})

#' @rdname setScale-method
#' @export
#' @keywords internal
getScale <- function(x) x@scale

#' @aliases setScale,PreprocessViews2-method
#' @rdname setScale-method
setReplaceMethod("setScale", "PreprocessViews2", function(x, value){
  x@scale <- value
  x
})

#' @aliases indexRanges,PreprocessViews2-method
#' @rdname indexRanges-method
setMethod("indexRanges", "PreprocessViews2", function(object) {
  rowRanges(object)$range_index
})

#' @rdname indexRanges-method
#' @aliases indexRanges<-,PreprocessViews2-method
setReplaceMethod("indexRanges", "PreprocessViews2", function(object, value) {
  rowRanges(object)$range_index <- value
  object
})

setMethod("paths", "PreprocessViews2", function(object) object@bamPaths)

setReplaceMethod("paths", "PreprocessViews2", function(object, value){
  object@bamPaths <- value
  object
})

#' @aliases rowRanges,PreprocessViews2,ANY-method
#' @rdname PreprocessViews2-class
setReplaceMethod("rowRanges", "PreprocessViews2", function(x, value){
  x@bamRanges <- value
  x
})

#' @aliases rowRanges,PreprocessViews2-method
#' @rdname PreprocessViews2-class
#' @param ... ignored
setMethod("rowRanges", "PreprocessViews2", function(x, ...) bamRanges(x))

#' Helper for creating filenames with .rds extension
#'
#' Intermediate files are stored in directories given by
#' \code{DataPaths}.  The names of the intermediate files are formed
#' by concatenating the \code{colnames} of the \code{BamViews}-derived
#' object with the extension \code{.rds}.
#'
#'
#' @examples
#'   library(Rsamtools)
#'   extdir <- system.file("extdata", package="Rsamtools", mustWork=TRUE)
#'   bamfile <- list.files(extdir, pattern="ex1.bam$", full.names=TRUE)
#'   bview <- BamViews(bamPaths=bamfile)
#'   rdsId(bview)
#'
#' 
#' @export
#' @param x a \code{BamViews}-derived object
rdsId <- function(x) {
  if(ncol(x) == 0) return(character())
  paste0(colnames(x), ".rds")
}

##--------------------------------------------------
##
## Coercion
##
##--------------------------------------------------

## setAs is exported in methods, so we do not export here

#' Coerce a \code{PreprocessViews2} object to a \code{RangedSummarizedExperiment}
#'
#' This method pulls the assay data from disk through the views object
#' interface, and then creates a \code{SummarizedExperiment} object
#' with an assay named 'copy'.
#'
#' @examples
#' pviews <- PreprocessViews2()
#' as(pviews, "RangedSummarizedExperiment")
#'
#' @return a \code{RangedSummarizedExperiment}
#' @param from character string ('PreprocessViews2')
#' @param to  character string  ('RangedSummarizedExperiment')
#' @rdname PreprocessViews2-coercion
#' @docType methods
#' @name setAs
#' @aliases coerce,PreprocessViews2,RangedSummarizedExperiment-method
setAs("PreprocessViews2", "RangedSummarizedExperiment", function(from, to){
  x <- assays(from)
  rr <- rowRanges(from)
  coldat <- bamSamples(from)
  SummarizedExperiment(assays=SimpleList(copy=x),
                       rowRanges=rr,
                       colData=coldat)
})


##--------------------------------------------------
##
## Descriptive stats 
##
##--------------------------------------------------

#' Compute the mean normalized read-depth for a set of intervals
#'
#' Calculates the mean normalized read-depth for a set of genomic
#' intervals in a \code{GRanges} object.
#'
#' @keywords internal
#' 
#' @export
#' @param pviews a \code{PreprocessViews2} object
#' @param gr  a \code{GRanges} object
granges_copynumber <- function(gr, pviews){
  pviews <- pviews[, 1]
  hits <- findOverlaps(gr, rowRanges(pviews))
  if(!any(duplicated(names(gr))) && !is.null(names(gr))){
    split_by <- factor(names(gr)[queryHits(hits)], levels=names(gr))
  } else {
    split_by <- queryHits(hits)
  }
  fc_context <- sapply(split(assays(pviews[subjectHits(hits), 1]), split_by), median)
  fc_context
}
