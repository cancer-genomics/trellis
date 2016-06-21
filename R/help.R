#' Functions to extract improper read pairs from bam files
#'
#'
#' @docType package
#' @name svalignments
#' @import methods
#' @import svclasses
#' @importMethodsFrom S4Vectors mcols mcols<- elementLengths 
#' @importFrom S4Vectors isTRUEorFALSE
#' @importFrom S4Vectors queryHits subjectHits
#' @importMethodsFrom GenomicAlignments readGAlignments first last 
#' @importFrom GenomicAlignments makeGAlignmentPairs findMateAlignment GAlignmentPairs GAlignments
#' @importFrom Rsamtools scanBamFlag bamPaths bamFlagAsBitMatrix 
#' @importMethodsFrom Rsamtools ScanBamParam 
#' @importMethodsFrom GenomeInfoDb seqnames seqlevels seqlevels<- seqinfo seqinfo<- seqlengths
#' @importMethodsFrom GenomeInfoDb seqlevelsStyle seqlevelsStyle<-
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom IRanges IRanges
#' @importMethodsFrom IRanges reduce width disjoin overlapsAny pintersect subsetByOverlaps unlist findOverlaps
#' @importFrom GenomicRanges GRanges GRangesList
#' @importMethodsFrom GenomicRanges granges strand strand<- width intersect 
#' @importMethodsFrom SummarizedExperiment rowRanges rowRanges<- 
#' @importFrom dplyr tbl_df
NULL
