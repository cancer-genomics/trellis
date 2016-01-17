#' Classes and methods for structural variation analyses
#'
#' svclass provides infrastructure for structural variant analyses of sequencing data.
#'
#' @docType package
#' @name svclasses
#' @import methods
#' @importMethodsFrom Rsamtools BamViews [ 
#' @importClassesFrom Rsamtools BamViews
#' @importFrom Rsamtools bamPaths bamRanges bamSamples
#' @importFrom S4Vectors SimpleList
#' @importMethodsFrom GenomeInfoDb seqnames seqlevels seqlevels<- seqinfo seqinfo<-
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importMethodsFrom GenomicRanges assays assays<- rowRanges rowRanges<-
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importMethodsFrom SummarizedExperiment SummarizedExperiment
NULL
