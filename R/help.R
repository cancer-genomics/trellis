#' Classes and methods for structural variation analyses
#'
#' svclass provides infrastructure for structural variant analyses of sequencing data.
#'
#' @docType package
#' @name svclasses
#' @import methods
#' @importFrom Rsamtools bamPaths bamRanges bamSamples
#' @importMethodsFrom GenomeInfoDb seqnames seqlevels
#' @importMethodsFrom Rsamtools BamViews [ 
#' @importClassesFrom Rsamtools BamViews
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importMethodsFrom GenomicRanges assays assays<- rowRanges rowRanges<-
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importMethodsFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors SimpleList
NULL
