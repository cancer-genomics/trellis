#' Classes and methods for structural variation analyses
#'
#' trellis provides infrastructure for structural variant analyses of sequencing data.
#'
#' @docType package
#' @name trellis
#' @import methods
#' @import BiocGenerics
#' @importMethodsFrom Rsamtools BamViews [
#' @importMethodsFrom Rsamtools ScanBamParam countBam pileup
#' @importClassesFrom Rsamtools BamViews
#' @importFrom Rsamtools bamPaths bamRanges bamSamples scanBamFlag PileupParam
#' @importFrom S4Vectors SimpleList mcols mcols<- queryHits subjectHits DataFrame
#' @importFrom S4Vectors queryHits subjectHits
#' @importMethodsFrom S4Vectors elementType elementNROWS
#' @importMethodsFrom GenomeInfoDb seqnames seqlevels seqlevels<- seqinfo seqinfo<- seqlengths
#' @importFrom GenomeInfoDb keepSeqlevels seqlevels seqinfo<- seqinfo
#' @importFrom IRanges IRanges LogicalList
#' @importMethodsFrom IRanges findOverlaps pintersect reduce overlapsAny subsetByOverlaps
#' @importFrom GenomicAlignments GAlignmentPairs GAlignments
#' @importClassesFrom GenomicAlignments GAlignmentPairs
#' @importMethodsFrom GenomicAlignments first last
#' @import GenomicRanges
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importMethodsFrom SummarizedExperiment SummarizedExperiment assays rowRanges
#' @importMethodsFrom SummarizedExperiment rowRanges<- colData colData<- assays<-
#' @importFrom matrixStats rowMedians
#' @importFrom graph graphNEL
#' @importMethodsFrom graph numNodes numEdges nodes addNode addEdge edges edgeNames
#' @importFrom DNAcopy CNA segment
#' @importFrom stats setNames predict acf loess median na.omit runif filter
#' @importFrom utils read.table data
#' @importFrom dplyr group_by
#' @importFrom magrittr '%>%'
NULL
