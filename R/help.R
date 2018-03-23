#' Somatic structural variant analysis
#'
#' trellis performs structural variant analyses of sequencing data.
#'
#' @docType package
#' @name trellis
#' @import methods
#' @import BiocGenerics
#' @importMethodsFrom Rsamtools BamViews [
#' @importMethodsFrom Rsamtools ScanBamParam countBam pileup scanBam
#' @importClassesFrom Rsamtools BamViews
#' @importFrom Rsamtools bamPaths bamRanges bamSamples scanBamFlag PileupParam scanBamWhat BamFile
#' @importFrom Rsamtools bamFlagAsBitMatrix bamMapqFilter
#' @importFrom S4Vectors SimpleList mcols mcols<- queryHits subjectHits DataFrame List
#' @importFrom S4Vectors queryHits subjectHits isTRUEorFALSE isRedundantHit isSelfHit
#' @importMethodsFrom S4Vectors elementType elementNROWS endoapply
#' @importMethodsFrom GenomeInfoDb seqnames seqlevels seqlevels<- seqinfo seqinfo<- seqlengths seqlengths<-
#' @importMethodsFrom GenomeInfoDb genome genome<- seqlevelsStyle<- seqlevelsStyle seqlevelsInUse
#' @importFrom GenomeInfoDb keepSeqlevels seqlevels seqinfo<- seqinfo
#' @importFrom IRanges IRanges LogicalList IntegerList %in%
#' @importMethodsFrom IRanges findOverlaps pintersect reduce overlapsAny subsetByOverlaps
#' @importMethodsFrom IRanges width disjoin pintersect unlist
#' @importFrom GenomicAlignments GAlignmentPairs GAlignments
#' @importFrom GenomicAlignments makeGAlignmentPairs findMateAlignment
#' @importClassesFrom GenomicAlignments GAlignmentPairs
#' @importMethodsFrom GenomicAlignments first last readGAlignments
#' @import GenomicRanges
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importMethodsFrom SummarizedExperiment SummarizedExperiment assays rowRanges
#' @importMethodsFrom SummarizedExperiment rowRanges<- colData colData<- assays<-
#' @importFrom matrixStats rowMedians rowMins
#' @importFrom graph graphNEL
#' @importMethodsFrom graph numNodes numEdges nodes addNode addEdge edges edgeNames
#' @importFrom DNAcopy CNA segment
#' @importFrom stats setNames predict acf loess median na.omit runif filter setNames kmeans
#' @importFrom utils read.table data read.delim head read.csv
#' @importFrom dplyr group_by top_n summarize mutate count n tbl_df
#' @importFrom magrittr '%>%'
#' @importFrom RBGL kCliques
#' @importFrom ggplot2 ggplot geom_rect ylab xlab scale_fill_manual aes
#' @importFrom ggplot2 theme scale_color_manual scale_x_continuous theme
#' @importFrom ggplot2 element_text element_blank facet_wrap
#' @importFrom scales trans_new pretty_breaks
#' @importFrom RColorBrewer brewer.pal
#' @importFrom AnnotationDbi select loadDb dbFileConnect dbFileDisconnect
#' @importMethodsFrom GenomicFeatures transcripts cdsBy exonsBy
#' @importFrom BSgenome getBSgenome
#' @importMethodsFrom Biostrings getSeq translate subseq
#' @importFrom Biostrings DNAStringSet
NULL
