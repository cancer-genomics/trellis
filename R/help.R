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
#' @importFrom S4Vectors queryHits subjectHits
#' @importMethodsFrom S4Vectors elementType elementNROWS
#' @importMethodsFrom GenomeInfoDb seqnames seqlevels seqlevels<- seqinfo seqinfo<- seqlengths
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom IRanges IRanges LogicalList
#' @importMethodsFrom IRanges findOverlaps pintersect reduce overlapsAny
#' @importFrom GenomicAlignments GAlignmentPairs GAlignments
#' @importClassesFrom GenomicAlignments GAlignmentPairs
#' @importMethodsFrom GenomicAlignments first last left right
#' @import GenomicRanges
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importMethodsFrom SummarizedExperiment SummarizedExperiment assays rowRanges 
#' @importMethodsFrom SummarizedExperiment rowRanges<- colData colData<- assays<-
#' @importFrom graph graphNEL 
#' @importMethodsFrom graph numNodes numEdges nodes addNode addEdge edges edgeNames
#' @importFrom stats median setNames
NULL
