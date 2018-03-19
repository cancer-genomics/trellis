<<<<<<< HEAD
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
#' @importMethodsFrom GenomicAlignments first last 
#' @import GenomicRanges
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importMethodsFrom SummarizedExperiment SummarizedExperiment assays rowRanges 
#' @importMethodsFrom SummarizedExperiment rowRanges<- colData colData<- assays<-
#' @importFrom graph graphNEL 
#' @importMethodsFrom graph numNodes numEdges nodes addNode addEdge edges edgeNames
#' @importFrom stats median setNames
=======
#' Preprocessing for whole genome structural variant analyses
#'
#' Preprocessing includes defining bin-level summaries of relative copy number,
#' modeling GC content, and excluding regions of the genome that are difficult
#' to interrogate for copy number alterations.
#'
#' @docType package
#' @name svpreprocess
#' @import methods
#' @import svclasses
#' @importFrom Rsamtools bamPaths bamRanges scanBamFlag
#' @importMethodsFrom Rsamtools ScanBamParam countBam
#' @importMethodsFrom IRanges overlapsAny findOverlaps
#' @importFrom S4Vectors subjectHits
#' @importMethodsFrom SummarizedExperiment rowRanges assays assays<-
#' @import GenomicRanges
#' @importFrom GenomeInfoDb keepSeqlevels seqlevels seqinfo<- seqinfo
#' @importFrom DNAcopy CNA segment
#' @importFrom stats setNames predict acf density loess mad median na.omit runif
>>>>>>> remote_svpreprocess/master
NULL
