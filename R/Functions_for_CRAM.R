#' Add seqinfo to dataframe 
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




#' Convert dataframe to GAlignment 
#'
#' @param dt 
#'
#' @return
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


#' properReadPairs
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
#' A broken-up version of properReadPairs_bedops that can be used in linux system. This
#' is the first part where we extract a set of GRanges that are near deletion regions. 
#' The second half of the properReadPairs_bedops function is properReadPairs_bedops2
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
#' Second half of the broken up properReadPairs_bedops function after we used bedmap 
#' of a bed file and the GRanges from prp_makeRthinned.
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


#' Define new functions for sv_deletions_bedops
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

#' First half of sv_deletion
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


#' Second half of sv_deletions
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