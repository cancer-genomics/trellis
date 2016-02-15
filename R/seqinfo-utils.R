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

#' @aliases chromosome,GAlignments-method
#' @rdname chromosome-methods
setMethod(chromosome, "GAlignments", function(object){
  as.character(seqnames(object))
})

#' @aliases chromosome,StructuralVariant-method
#' @rdname chromosome-methods
setMethod("chromosome", "StructuralVariant", function(object) chromosome(variant(object)))

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

#' @param x a \code{PreprocessViews2} object
#' @aliases seqlevels,PreprocessViews2-method
#' @rdname PreprocessViews2-class
setMethod("seqlevels", "PreprocessViews2", function(x) seqlevels(rowRanges(x)))

#' @rdname RearrangementList-class
#' @aliases seqlevels,RearrangementList-method
setMethod("seqlevels", "RearrangementList", function(x){
  seqlevels(linkedBins(x[[1]]))
})


#' @aliases seqinfo,PreprocessViews2-method
#' @rdname PreprocessViews2-class
setMethod("seqinfo", "PreprocessViews2", function(x) seqinfo(rowRanges(x)))

#' @param value an object of class \code{seqinfo}
#' @aliases seqinfo<-,PreprocessViews2,seqinfo-method
#' @rdname PreprocessViews2-class
setReplaceMethod("seqinfo", "PreprocessViews2",
                 function(x, value){
                   rowRanges(x) <- keepSeqlevels(rowRanges(x), seqlevels(value))
                   ##seqlevels(rowRanges(x), force=TRUE) <- seqlevels(value)
                   seqinfo(rowRanges(x)) <- value
                   x
                 })

#' @aliases seqlengths,StructuralVariant-method
#' @rdname StructuralVariant-class
setMethod("seqlengths", "StructuralVariant", function(x) seqlengths(variant(x)))

#' @aliases seqinfo,StructuralVariant-method
#' @rdname StructuralVariant-class
setMethod("seqinfo", "StructuralVariant", function(x) seqinfo(variant(x)))
