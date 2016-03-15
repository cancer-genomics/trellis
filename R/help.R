#' Identify copy number variants from log2-transformed and GC-adjusted read counts.
#'
#' Segments bin-level estimates of copy number.  Use both read pairs
#' and read depth to identify focal deletions and amplicons. For
#' amplicons, we use read pairs to identify grouped amplicons.
#'
#' @docType package
#' @name svcnvs
#' @import methods
#' @import svclasses
#' @import S4Vectors
#' @importFrom BiocGenerics combine
#' @importFrom Rsamtools bamPaths
#' @importFrom DNAcopy CNA segment
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importMethodsFrom GenomeInfoDb seqlevelsStyle<- seqlevelsStyle seqlevels seqlevels<- seqlengths
#' @importMethodsFrom GenomeInfoDb seqlevels seqnames seqinfo
#' @importFrom GenomicAlignments GAlignmentPairs
#' @importMethodsFrom GenomicAlignments first last
#' @importClassesFrom GenomicAlignments GAlignmentPairs
#' @importFrom IRanges IRanges
#' @import GenomicRanges
#' @importFrom svalignments properReadPairs readPairsNearVariant get_readpairs get_improper_readpairs totalWidth
#' @importFrom graph graphNEL
#' @importMethodsFrom graph numNodes numEdges edges addEdge addNode nodes
#' @importFrom RBGL kCliques
NULL
