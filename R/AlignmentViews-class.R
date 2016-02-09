#' @include Preprocess-class.R
NULL

#' Container for storing views to extracted alignments
#'
#' Improper read pairs are extracted from BAM files to augment CNV-
#' and rearrangement-analyses
#' 
#' @slot improperAlignmentPaths character-vector of file paths 
#' @export
#' @rdname AlignmentViews2-class
setClass("AlignmentViews2", contains="PreprocessViews2",
         representation(improperAlignmentPaths="character"))


setValidity("AlignmentViews2", function(object){
  msg <- TRUE
  if(length(improperPaths(object)) != ncol(object)){
    msg <- "length of improperPaths must be the same as ncol"
    return(msg)
  }
  msg
})


#' Constructor for AlignmentViews2 objects
#'
#' AlignmentViews2 objects inherit from \code{BamViews}.  In addition
#' to containing file paths to one or more bam files, AlignmentViews2
#' objects contain file paths to a subset of improperly paired reads
#' that were / will be parsed from the bams.  There will be one parsed
#' file per bam file.
#'
#' @examples
#'   library(Rsamtools)
#'   require(TestBams)
#'   extdata <- system.file("extdata", package="TestBams")
#'   bam.file <- list.files(extdata, pattern="\\.bam$", full.name=TRUE)
#'   bv <- BamViews(bam.file)
#'   dp <- DataPaths(tempdir())
#'   aviews <- AlignmentViews2(bv, dp)
#'   validObject(aviews[, 1])
#'   aviews
#'   improperPaths(aviews)
#' @export
#' @param bviews A \code{BamViews} object
#' @param dirs A \code{DataPaths} object
#' @rdname AlignmentViews2-class
AlignmentViews2 <- function(bviews, dirs){
  path <- dirs[["alignments/0improper"]]
  improper_paths <- file.path(path, rdsId(bviews))
  aviews <- as(bviews, "AlignmentViews2")
  indexRanges(aviews) <- seq_len(nrow(aviews))
  improperPaths(aviews) <- improper_paths
  aviews
}

#' Accessor for file names containing improper read pair alignments
#'
#' @return a character-vector of file paths
#' @rdname AlignmentViews2-class
#' @export
#' @param object a \code{AlignmentViews2} object
setGeneric("improperPaths", function(object) standardGeneric("improperPaths"))

setGeneric("improperPaths<-", function(object, value) standardGeneric("improperPaths<-"))


#' @aliases improperPaths,AlignmentViews2-method
#' @rdname AlignmentViews2-class
setMethod("improperPaths", "AlignmentViews2", function(object)
  object@improperAlignmentPaths)

setReplaceMethod("improperPaths", "AlignmentViews2",
                 function(object, value){
                   object@improperAlignmentPaths <- value
                   object
                 })

##--------------------------------------------------
##
## Subsetting
##
##--------------------------------------------------

#' @rdname AlignmentViews2-class
#' @aliases [,AlignmentViews2,missing,numeric-method
#' @param x a \code{AlignmentViews2} object
#' @param i index for \code{rowRanges} elements 
#' @param j index for samples
#' @param ... ignored
#' @param drop ignored
setMethod("[", c("AlignmentViews2", "missing", "numeric"), function(x, i, j, ..., drop=FALSE){
  improperPaths(x) <- improperPaths(x)[j]
  x <- callNextMethod(x, i, j, ..., drop)
  x
})

#' @rdname AlignmentViews2-class
#' @aliases [,AlignmentViews2,missing,logical-method
setMethod("[", c("AlignmentViews2", "missing", "logical"), function(x, i, j, ..., drop=FALSE){
  improperPaths(x) <- improperPaths(x)[j]
  x <- callNextMethod(x, i, j, ..., drop)
  x
})

#' @rdname AlignmentViews2-class
#' @aliases [,AlignmentViews2,missing,character-method
setMethod("[", c("AlignmentViews2", "missing", "character"), function(x, i, j, ..., drop=FALSE){
  k <- match(j, colnames(x))
  improperPaths(x) <- improperPaths(x)[k]
  x <- callNextMethod(x, i, j, ..., drop)
  x
})

