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
#' @import S4Vectors
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
#' @importFrom stats setNames predict acf loess median runif setNames kmeans
#' @importFrom utils read.table data read.delim read.csv
#' @importFrom tibble as.tibble
#' @importFrom tidyr unite
#' @importFrom dplyr group_by top_n summarize mutate count n tbl_df filter left_join ungroup
#' @importFrom purrr map_lgl
#' @importFrom magrittr '%>%' set_names
#' @importFrom RBGL kCliques
#' @importFrom ggplot2 ggplot geom_rect ylab xlab scale_fill_manual aes
#' @importFrom ggplot2 theme scale_color_manual scale_x_continuous theme coord_cartesian guides geom_vline guide_legend ggtitle scale_x_reverse
#' @importFrom ggplot2 element_text element_blank facet_wrap ggplotGrob
#' @importFrom scales trans_new pretty_breaks
#' @importFrom RColorBrewer brewer.pal
#' @importFrom AnnotationDbi select loadDb dbFileConnect dbFileDisconnect
#' @importMethodsFrom GenomicFeatures transcripts cdsBy exonsBy
#' @importFrom BSgenome getBSgenome
#' @importMethodsFrom Biostrings getSeq translate subseq
#' @importFrom Biostrings DNAStringSet
#' @importFrom gridExtra grid.arrange
#' @importFrom grid unit
NULL
