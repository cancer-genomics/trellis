#' Add seqinfo to dataframe (for CRAM)
#'
#' @param data 
#'
#' @return data with added seqinfo 
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




#' Convert dataframe to GAlignment (for CRAM)
#'
#' @param dt 
#'
#' @return GAlignment object that is ready to be paired with GAlignmentPairs
#' @export
dt2ga2 <- function(dt){
  # sort file by chromosome in alpha-numeric order
  #dt <- dt[gtools::mixedorder(dt$chr), ]
  
  # seqlengths()
  chrom.info <- as.data.frame(seqinfo(Hsapiens))[1:23,]
  chrom <- levels(Rle(factor(dt$chr)))
  seq.info <- chrom.info[chrom, 1]
  seqlen <- setNames(seq.info, chrom)
  
  # if seq is present in data frame 
  if ("seq" %in% colnames(dt)){
    seq = as.character(dt$seq)
  }
  else {seq = NA}
  
  ga <- GAlignments(seqnames = Rle(factor(dt$chr)), 
                    pos = as.integer(dt$start+1), 
                    cigar = as.character(dt$cigar), 
                    strand = GenomicRanges::strand(dt$strand), 
                    names = as.character(dt$id), 
                    seqlengths = seqlen, 
                    flag = as.integer(dt$flag), 
                    mrnm = as.factor(dt$mrnm),
                    mpos = as.integer(dt$mpos),
                    mapq = as.numeric(dt$score), 
                    seq = seq, 
                    qname = as.character(dt$id))
  ga
} # +1 in the starting position since bed files are 
#' @examples




## Functions for Deletion


#' Proper read pairs near variant (properReadPairs for CRAM)
#'
#' @param prp4del.ga 
#' @param del.gr 
#' @param dp 
#'
#' @return A GAlignmentPairs object that are proper reads near the deletion region
#' @export
properReadPairs_bedops <- function(prp4del.ga, del.gr, dp){
  #Extract GRanges to be used as query 
  thinned <- del.gr %>% 
    thinReadPairQuery() %>%
    filter(start > 0 & end > 0) # some fragments had negative values in start and end
  seqlevelsStyle(thinned) <- bamSeqLevelsStyle(dp)
  rthinned <- IRanges::reduce(thinned) # remove overlaps, merge overlaps into one large fragment
  
  ga <- subsetByOverlaps(prp4del.ga, rthinned)
  galp <- makeGAlignmentPairs(ga, use.names=TRUE, use.mcols=TRUE, strandMode=1)
  is_valid <- validPairForDeletion(galp)
  proper_rp.bedops <- galp[is_valid]
  proper_rp.bedops
}

# break-up of properReadPairs_bedops
# first part
#' First part of properReadPairs_bedops 
#' 
#' This is the first part of the broken-up version of properReadPairs_bedops that can 
#' be used in linux system. This is the first part where we extract a set of GRanges 
#' that are near deletion regions. The second half of the properReadPairs_bedops function is properReadPairs_bedops2
#'
#' @param del.gr 
#' @param dp 
#'
#' @return GRanges that are near deletion regions
#' @export
prp_makeRthinned <- function(del.gr, dp){
  thinned <- del.gr %>% 
    thinReadPairQuery() %>%
    filter(start > 0 & end > 0) # some fragments had negative values in start and end
  seqlevelsStyle(thinned) <- bamSeqLevelsStyle(dp)
  rthinned <- IRanges::reduce(thinned)
  rthinned
}


# second part: after bedmap in shell file, import GAlignment
#' Second part of properReadPairs_bedops 
#' 
#' Used after we used bedmap of a bed file and the GRanges from prp_makeRthinned 
#'
#' @param ga 
#'
#' @return A set of proper read pairs near deletion region
#' @export
properReadPairs_bedops2 <- function(ga){
  galp <- makeGAlignmentPairs(ga, use.names=TRUE, use.mcols=TRUE, strandMode=1)
  is_valid <- validPairForDeletion(galp)
  proper_rp.bedops <- galp[is_valid]
  proper_rp.bedops
}


# sv_deletions
readPairsNearVariant_bedops <- function(object, query){
  ga <- subsetByOverlaps(object, query)
  galp <- makeGAlignmentPairs(ga, use.names=TRUE, use.mcols=TRUE, strandMode=1)
  seqlevelsStyle(galp) <- seqlevelsStyle(query)
  is_valid <- validPairForDeletion(galp)
  proper_rp <- galp[is_valid]
  proper_rp
}

allProperReadPairs_bedops <- function(sv, param, object, zoom.out=1){
  cnv <- variant(sv)
  g <- expandGRanges(cnv, width(cnv)*zoom.out)
  reduce <- IRanges::reduce
  query <- reduce(disjoin(g))
  proper_rp <- readPairsNearVariant_bedops(object, query)
  sv@proper <- proper_rp
  indexProper(sv) <- initializeProperIndex3(sv, zoom.out=1)
  sv
}

finalize_deletions_bedops <- function(sv, preprocess, gr_filters,
                                      param=DeletionParam(), object){
  if(length(sv) == 0) return(sv)
  if(missing(gr_filters)){
    gr_filters <- genomeFilters(preprocess$genome)
  }
  bins <- preprocess$bins
  sv <- SVFilters(sv, gr_filters, bins, param=param)
  if (length(sv) == 0) return(sv)
  sv <- groupSVs(sv)
  
  sv <- allProperReadPairs_bedops(sv, param, object, zoom.out=1)
  
  if(length(sv@proper) > 25e3){
    proper(sv) <- sv@proper[sample(seq_along(sv@proper), 25e3)]
    indexProper(sv) <- initializeProperIndex3(sv, zoom.out=1)
  }
  sv <- rename(sort(sv))
  sv
}


#' Define new functions for sv_deletions (for CRAM)
#'
#' @param preprocess 
#' @param gr_filters 
#' @param param 
#' @param object 
#'
#' @return A set of candidate structural variants
#' @export
sv_deletions_bedops <- function(preprocess, gr_filters, param=DeletionParam(), object){
  if(missing(gr_filters)){
    gr_filters <- genomeFilters(preprocess$genome)
  }
  
  segs2 <- preprocess$segments
  preprocess$segments <- segs2[segs2$seg.mean < hemizygousThr(param)]
  sv <- deletion_call(preprocess, gr_filters, param)
  calls(sv) <- rpSupportedDeletions(sv, param=param, bins=preprocess$bins)
  
  if(param@remove_hemizygous){
    sv <- removeHemizygous(sv)
  }
  # Extracting improper read pairs with a high MAPQ
  improper_rp <- preprocess$read_pairs[["improper"]]
  ## avoid namespace issues with dplyr
  first <- GenomicAlignments::first
  last <- GenomicAlignments::last
  mapq <- mcols(first(improper_rp))$mapq > 30 &
    mcols(last(improper_rp))$mapq > 30
  improper_rp <- improper_rp[mapq]
  
  if(length(variant(sv)) > 0){
    # Revise deletion boundaries using improper read pairs
    # and proper read pairs 
    sv <- reviseEachJunction(sv=sv, bins=preprocess$bins,
                             improper_rp=improper_rp,
                             param=param)
    if(param@remove_hemizygous) {
      sv <- removeHemizygous(sv)
    }
  }
  
  # stop here since finalize_deletions() requires bam file
  
  if(length(variant(sv)) > 0){
    sv <- revise(sv, bins=preprocess$bins, param=param)
    sv <- finalize_deletions_bedops(sv=sv, preprocess,
                                    gr_filters=gr_filters,
                                    param=param, 
                                    object=object)
  }
  sv
}



# break up of sv_deletions -  make query file 
# first half of allProperReadPairs_bedops for making query
allProperReadPairs_bedops1 <- function(sv, zoom.out=1){
  cnv <- variant(sv)
  g <- expandGRanges(cnv, width(cnv)*zoom.out)
  reduce <- IRanges::reduce
  query <- reduce(disjoin(g))
  query
}


# first half of finalize_deletions_bedops for making query
finalize_deletions_bedops1 <- function(sv, preprocess, gr_filters,
                                       param=DeletionParam()){
  if(length(sv) == 0) return(sv)
  if(missing(gr_filters)){
    gr_filters <- genomeFilters(preprocess$genome)
  }
  bins <- preprocess$bins
  sv <- SVFilters(sv, gr_filters, bins, param=param)
  if (length(sv) == 0) return(sv)
  sv <- groupSVs(sv)
  query <- allProperReadPairs_bedops1(sv, zoom.out=1)
  query
}

#' First part of sv_deletions_bedops
#'
#' @param preprocess 
#' @param gr_filters 
#' @param param 
#'
#' @return A GRanges object that can be used for bedmap in command-line programs
#' @export
sv_del_makeQuery <- function(preprocess, gr_filters, param=DeletionParam()){
  # Part 1: first half of sv_deletions_bedops
  if(missing(gr_filters)){
    gr_filters <- genomeFilters(preprocess$genome)
  }
  
  segs2 <- preprocess$segments
  preprocess$segments <- segs2[segs2$seg.mean < hemizygousThr(param)]
  sv <- deletion_call(preprocess, gr_filters, param)
  calls(sv) <- rpSupportedDeletions(sv, param=param, bins=preprocess$bins)
  
  if(param@remove_hemizygous){
    sv <- removeHemizygous(sv)
  }
  # Extracting improper read pairs with a high MAPQ
  improper_rp <- preprocess$read_pairs[["improper"]]
  ## avoid namespace issues with dplyr
  first <- GenomicAlignments::first
  last <- GenomicAlignments::last
  mapq <- mcols(first(improper_rp))$mapq > 30 &
    mcols(last(improper_rp))$mapq > 30
  improper_rp <- improper_rp[mapq]
  
  if(length(variant(sv)) > 0){
    # Revise deletion boundaries using improper read pairs
    # and proper read pairs 
    sv <- reviseEachJunction(sv=sv, bins=preprocess$bins,
                             improper_rp=improper_rp,
                             param=param)
    if(param@remove_hemizygous) {
      sv <- removeHemizygous(sv)
    }
  }
  
  # stop here since finalize_deletions() requires bam file
  
  if(length(variant(sv)) > 0){
    sv <- revise(sv, bins=preprocess$bins, param=param)
    # Part 2: fist half of finalize_deletions_bedops
    query <- finalize_deletions_bedops1(sv=sv, preprocess,
                                        gr_filters=gr_filters,
                                        param=param)
  }
  query
}


#second half of sv_deletions_bedops after bedmap
readPairsNearVariant_bedops2 <- function(object, query){
  #ga <- subsetByOverlaps(object, query)
  galp <- makeGAlignmentPairs(object, use.names=TRUE, use.mcols=TRUE, strandMode=1)
  seqlevelsStyle(galp) <- seqlevelsStyle(query)
  is_valid <- validPairForDeletion(galp)
  proper_rp <- galp[is_valid]
  proper_rp
}

allProperReadPairs_bedops2 <- function(sv, param, object, zoom.out=1){
  cnv <- variant(sv)
  g <- expandGRanges(cnv, width(cnv)*zoom.out)
  reduce <- IRanges::reduce
  query <- reduce(disjoin(g))
  proper_rp <- readPairsNearVariant_bedops2(object, query)
  sv@proper <- proper_rp
  indexProper(sv) <- initializeProperIndex3(sv, zoom.out=1)
  sv
}

finalize_deletions_bedops2 <- function(sv, preprocess, gr_filters,
                                       param=DeletionParam(), object){
  if(length(sv) == 0) return(sv)
  if(missing(gr_filters)){
    gr_filters <- genomeFilters(preprocess$genome)
  }
  bins <- preprocess$bins
  sv <- SVFilters(sv, gr_filters, bins, param=param)
  if (length(sv) == 0) return(sv)
  sv <- groupSVs(sv)
  
  sv <- allProperReadPairs_bedops2(sv, param, object, zoom.out=1)
  
  if(length(sv@proper) > 25e3){
    proper(sv) <- sv@proper[sample(seq_along(sv@proper), 25e3)]
    indexProper(sv) <- initializeProperIndex3(sv, zoom.out=1)
  }
  sv <- rename(sort(sv))
  sv
}


#' Second part of sv_deletions_bedops
#'
#' @param preprocess 
#' @param gr_filters 
#' @param param 
#' @param object 
#'
#' @return A set of candidate structural variants.  
#' @export
sv_deletions_bedops2 <- function(preprocess, gr_filters, param=DeletionParam(), object){
  if(missing(gr_filters)){
    gr_filters <- genomeFilters(preprocess$genome)
  }
  
  segs2 <- preprocess$segments
  preprocess$segments <- segs2[segs2$seg.mean < hemizygousThr(param)]
  sv <- deletion_call(preprocess, gr_filters, param)
  calls(sv) <- rpSupportedDeletions(sv, param=param, bins=preprocess$bins)
  
  if(param@remove_hemizygous){
    sv <- removeHemizygous(sv)
  }
  # Extracting improper read pairs with a high MAPQ
  improper_rp <- preprocess$read_pairs[["improper"]]
  ## avoid namespace issues with dplyr
  first <- GenomicAlignments::first
  last <- GenomicAlignments::last
  mapq <- mcols(first(improper_rp))$mapq > 30 &
    mcols(last(improper_rp))$mapq > 30
  improper_rp <- improper_rp[mapq]
  
  if(length(variant(sv)) > 0){
    # Revise deletion boundaries using improper read pairs
    # and proper read pairs 
    sv <- reviseEachJunction(sv=sv, bins=preprocess$bins,
                             improper_rp=improper_rp,
                             param=param)
    if(param@remove_hemizygous) {
      sv <- removeHemizygous(sv)
    }
  }
  
  # stop here since finalize_deletions() requires bam file
  
  if(length(variant(sv)) > 0){
    sv <- revise(sv, bins=preprocess$bins, param=param)
    sv <- finalize_deletions_bedops2(sv=sv, preprocess,
                                     gr_filters=gr_filters,
                                     param=param, 
                                     object=object)
  }
  sv
}



## Amplification 

# Define functions for sv_amplicons2_bedops
get_readpairs2_bedops <- function(g=queryRanges(ag), bed_file){
  #g <- queryRanges(ag.bedops)
  ga <- subsetByOverlaps(bed_file, g)
  galp <- makeGAlignmentPairs(ga, use.names=TRUE, use.mcols=TRUE, strandMode=1)
  validR1 <- overlapsAny(GenomicAlignments::first(galp), g)
  validR2 <- overlapsAny(GenomicAlignments::last(galp), g)
  proper_rp <- galp[validR1 & validR2]
  proper_rp
}

add_amplicons_bedops <- function(ag, bed_file, params){
  proper_rp <- get_readpairs2_bedops(g=queryRanges(ag), bed_file)
  ag <- addFocalDupsFlankingAmplicon(ag, proper_rp, params)
  queryRanges(ag) <- focalAmpliconDupRanges(ag, params)
  ag
}


#' sv_amplicons2 for CRAM
#'
#' @param preprocess 
#' @param amplicon_filters 
#' @param params 
#' @param bed_file 
#'
#' @return A set of amplicons
#' @export
sv_amplicons2_bedops <- function(preprocess, amplicon_filters, params=ampliconParams(), bed_file=amp.ga){
  if(missing(amplicon_filters)){
    amplicon_filters <- ampliconFilters(preprocess$genome)
  }
  preprocess$segments <- amplified_segments(preprocess$segments, params)
  ag <- initialize_graph2(preprocess, amplicon_filters, params)
  if (length(ag) == 0) return(ag)
  
  ag <- add_amplicons_bedops(ag, bed_file, params)
  
  improper_rp <- preprocess$read_pairs[["improper"]]
  ag <- link_amplicons(ag, improper_rp, params)
  tx <- loadTx(preprocess$genome)
  ag <- annotate_amplicons(ag, tx)
  ag
}


#sv_amplicons2_bedops first half

# output amplicon query
#' First part of sv_amplicons2_bedops
#'
#' @param preprocess 
#' @param amplicon_filters 
#' @param params 
#'
#' @return A set of GRanges near amplicons
#' @export
sv_amp_makeQuery <- function(preprocess, amplicon_filters, params=ampliconParams()){
  if(missing(amplicon_filters)){
    amplicon_filters <- ampliconFilters(preprocess$genome)
  }
  preprocess$segments <- amplified_segments(preprocess$segments, params)
  ag <- initialize_graph2(preprocess, amplicon_filters, params)
  if (length(ag) == 0) return(ag)
  
  g=queryRanges(ag)
  g
}


#sv_amplicons2_bedops second half
# Define functions for sv_amplicons2_bedops
get_readpairs2_bedops2 <- function(g=queryRanges(ag), bed_file){
  #g <- queryRanges(ag.bedops)
  #ga <- subsetByOverlaps(bed_file, g)
  galp <- makeGAlignmentPairs(bed_file, use.names=TRUE, use.mcols=TRUE, strandMode=1)
  validR1 <- overlapsAny(GenomicAlignments::first(galp), g)
  validR2 <- overlapsAny(GenomicAlignments::last(galp), g)
  proper_rp <- galp[validR1 & validR2]
  proper_rp
}

add_amplicons_bedops2 <- function(ag, bed_file, params){
  proper_rp <- get_readpairs2_bedops2(g=queryRanges(ag), bed_file)
  ag <- addFocalDupsFlankingAmplicon(ag, proper_rp, params)
  queryRanges(ag) <- focalAmpliconDupRanges(ag, params)
  ag
}


#' Second part of sv_amplicons2_bedops
#'
#' @param preprocess 
#' @param amplicon_filters 
#' @param params 
#' @param bed_file 
#'
#' @return A set of amplicons
#' @export
sv_amplicons2_bedops2 <- function(preprocess, amplicon_filters, params=ampliconParams(), bed_file){
  if(missing(amplicon_filters)){
    amplicon_filters <- ampliconFilters(preprocess$genome)
  }
  preprocess$segments <- amplified_segments(preprocess$segments, params)
  ag <- initialize_graph2(preprocess, amplicon_filters, params)
  if (length(ag) == 0) return(ag)
  
  ag <- add_amplicons_bedops2(ag, bed_file, params)
  
  improper_rp <- preprocess$read_pairs[["improper"]]
  ag <- link_amplicons(ag, improper_rp, params)
  tx <- loadTx(preprocess$genome)
  ag <- annotate_amplicons(ag, tx)
  ag
}




#Rearrangement
#BLAT mapped-mapped

# for mapped-mapped, BLAT realignment
# define functions to convert GAlignment to tibble for first and last reads
ga2tibble_first <- function(first, len){
  first.tib <- as_tibble(first)
  rname <- as.character(seqnames(first))
  read <- rep("R1", len)
  id <- rep(paste0(sample, ".bam"), len)
  rearrangement.id <- names(first)
  first.tib <- first.tib %>% dplyr::select(-rpid) %>% mutate(rname, read, id, rearrangement.id)
  first.tib
}

ga2tibble_last <- function(last, len){
  last.tib <- as_tibble(last)
  rname <- as.character(seqnames(last))
  read <- rep("R2", len)
  id <- rep(paste0(sample, ".bam"), len)
  rearrangement.id <- names(last)
  last.tib <- last.tib %>% dplyr::select(-rpid) %>% mutate(rname, read, id, rearrangement.id)
  last.tib
}

#' getSequenceOfReads for CRAM
#'
#' @param rlist.bedops 
#' @param MAX 
#' @param sample 
#'
#' @return Tagged sequences from BLAT
#' @export
getSequenceOfReads_bedops <- function(rlist.bedops, MAX=25, sample=sample){
  # initialize an empty tibble
  #colnames(tags)
  tags_colnames <- c("seqnames", "strand", "cigar", "qwidth", "start", "end", "width", "njunc", "qname", "rname", "seq", "flag", "mrnm", "mpos", "mapq", "read", "id", "rearrangement.id")
  tags_bedops <- as_tibble(matrix(nrow = 0, ncol = length(tags_colnames)), .name_repair = ~ tags_colnames)
  
  for (i in 1:length(rlist.bedops)) {
    imp <- improper(rlist.bedops[[i]])
    len <- length(imp)
    
    if (len > MAX) {
      rown <- sample(1:len, MAX)
      imp_sub <- imp[rown, ]
      first <- GenomicAlignments::first(imp_sub)
      last <- GenomicAlignments::last(imp_sub)
      len <- MAX
    }
    else {
      first <- GenomicAlignments::first(imp)
      last <- GenomicAlignments::last(imp)
    }
    
    # compile tibble
    first_tib <- ga2tibble_first(first, len)
    last_tib <- ga2tibble_last(last, len)
    tib <- rbind(first_tib, last_tib)
    
    # append to tags tibble 
    tags_bedops <- rbind(tags_bedops, tib)
  }
  tags_bedops
}



#BLAT mapped-unmapped

# function to import R1 and R2 as granges object
import_mum_bedops <- function(path, sample=sample, file){
  mum.df <- data.table::fread(paste0(path, "/", sample, "_", file))
  colnames(mum.df) <- c("chr", "start", "end", "snms", "strand", "seq")
  # replace all strand values with *
  mum.df$strand <- rep("*", dim(mum.df)[1])
  # convert data.table to GRanges
  # bed file start at +1 compared to non-bed file -- use starts.in.df.are.0based = TRUE to be the same as vignette
  mum.gr <- GenomicRanges::makeGRangesFromDataFrame(mum.df, starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
  names(mum.gr) <- mum.df$snms
  mum.gr
}

#' unmapped_read for CRAM
#'
#' @param query.bedops 
#' @param maxgap 
#' @param path 
#' @param sample 
#'
#' @return Unmapped reads for BLAT
#' @export
unmapped_read_bedops <- function(query.bedops, maxgap=500, path, sample=sample){
  # import R1 and R2 as granges object
  mumR1.name <- "mum_R1.bed"
  mumR2.name <- "mum_R2.bed"
  mumR1.gr <- import_mum_bedops(path, sample = sample, file = mumR1.name)
  mumR2.gr <- import_mum_bedops(path, sample = sample, file = mumR2.name)
  
  mate_gr <- subsetByOverlaps(mumR1.gr, query.bedops+maxgap)
  mate_gr$read <- rep("R1", length(mate_gr))
  mate_gr2 <- subsetByOverlaps(mumR2.gr, query.bedops+maxgap)
  mate_gr2$read <- rep("R2", length(mate_gr2))
  unmapped.bedops <- c(mate_gr, mate_gr2)
  unmapped.bedops
}



