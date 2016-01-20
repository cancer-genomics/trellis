#' Functions to extract improper read pairs from bam files
#'
#'
#' @docType package
#' @name svalignments
#' @import methods
#' @import svclasses
#' @importMethodsFrom S4Vectors mcols mcols<-
#' @importFrom S4Vectors isTRUEorFALSE
#' @importMethodsFrom GenomicAlignments readGAlignments first last
#' @importFrom GenomicAlignments makeGAlignmentPairs findMateAlignment GAlignmentPairs
#' @importFrom Rsamtools scanBamFlag bamPaths bamFlagAsBitMatrix 
#' @importMethodsFrom Rsamtools ScanBamParam 
#' @importMethodsFrom GenomeInfoDb seqnames seqlevels seqlevels<- seqinfo seqinfo<- seqlengths seqlevelsStyle seqlevelsStyle<-
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom IRanges IRanges
#' @importMethodsFrom IRanges reduce width disjoin
#' @importFrom GenomicRanges GRanges
#' @importMethodsFrom GenomicRanges granges strand
NULL
