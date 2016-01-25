#' Classes and methods for structural variation analyses
#'
#' svclass provides infrastructure for structural variant analyses of sequencing data.
#'
#' @docType package
#' @name svclasses
#' @import methods
#' @import BiocGenerics
#' @importMethodsFrom Rsamtools BamViews [ 
#' @importClassesFrom Rsamtools BamViews
#' @importFrom Rsamtools bamPaths bamRanges bamSamples
#' @importFrom S4Vectors SimpleList mcols mcols<- queryHits subjectHits DataFrame
#' @importMethodsFrom S4Vectors elementType queryHits subjectHits elementLengths
#' @importMethodsFrom GenomeInfoDb seqnames seqlevels seqlevels<- seqinfo seqinfo<- seqlengths
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom IRanges IRanges
#' @importMethodsFrom IRanges findOverlaps pintersect reduce
#' @importFrom GenomicAlignments GAlignmentPairs GAlignments
#' @importClassesFrom GenomicAlignments GAlignmentPairs
#' @import GenomicRanges
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importMethodsFrom SummarizedExperiment SummarizedExperiment
NULL
