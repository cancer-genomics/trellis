#' Extract autosomal (human) sequence levels
#' 
#' @export
#' @param object a \code{GRanges} object
#' @docType methods
#' @rdname chromosome-methods
setGeneric("autosomeNames", function(object) standardGeneric("autosomeNames"))

#' Return seqnames of a \code{GRanges} object as a character vector
#' 
#' @export
#' @docType methods
#' @rdname chromosome-methods
setGeneric("chromosome", function(object) standardGeneric("chromosome"))

#' @aliases chromosome,GRanges-method
#' @rdname chromosome-methods
setMethod(chromosome, "GRanges", function(object) as.character(seqnames(object)))

#' @aliases chromosome,RangedSummarizedExperiment-method
#' @rdname chromosome-methods
setMethod(chromosome, "RangedSummarizedExperiment", function(object){
  chromosome(rowRanges(object))
})

#' @aliases autosomeNames,BamViews-method
#' @rdname chromosome-methods
setMethod(autosomeNames, "BamViews", function(object){
  autosomeNames(bamRanges(object))
})

#' @aliases autosomeNames,GRanges-method
#' @rdname chromosome-methods
setMethod(autosomeNames, "GRanges", function(object){
    sl <- seqlevels(object)
    sl[!sl %in% c("chrX", "chrY", "chrM", "X", "Y", "M")]
})

##--------------------------------------------------
##
## PreprocessViews2 methods
##
##--------------------------------------------------

setMethod("seqlevels", "PreprocessViews2", function(x) seqlevels(rowRanges(x)))

setMethod("seqinfo", "PreprocessViews2", function(x) seqinfo(rowRanges(x)))

setReplaceMethod("seqinfo", "PreprocessViews2",
                 function(x, value){
                   seqlevels(rowRanges(x), force=TRUE) <- seqlevels(value)
                   seqinfo(rowRanges(x)) <- value
                   x
                 })
