#' Add seqinfo to bed files (for cram file version)
#'
#' @param data 
#'
#' @return data
#' @export
add.info <- function(data){
  chrom.info <- as.data.frame(seqinfo(Hsapiens))[1:23,]
  chrom <- levels(seqnames(data))
  #seq.len <- chrom.info[chrom, 1]
  #seqlen <- setNames(seq.info, chrom)
  seqlengths(data) <- chrom.info[chrom,1]
  isCircular(data) <- chrom.info[chrom,2]
  genome(data) <- chrom.info[chrom,3]
  data
}

#'
#' @examples
#'   library(svfilters.hg19)
#'   data(bins1kb)
#'   library(GenomeInfoDb)
#'   library(DNAcopy)
#'   fc.bowtie <- data.table::fread(here("cgov44t_comparison", "bowtie", "cgov44t_revised_fragcount.bed"))
#'   colnames(fc.bowtie) <- c("chr", "start", "end", "gc", "map", "cnt")
#'   bins1kb.bowtie <- GenomicRanges::makeGRangesFromDataFrame(fc.bowtie, starts.in.df.are.0based = FALSE, keep.extra.columns = TRUE) #starts.in.df.are.0based = TRUE makes ranges start from 01 not 02
#'   bins.bedops.bowtie <- add.info(bins1kb.bowtie)