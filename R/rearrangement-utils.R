#' @include AllGenerics.R
NULL

##--------------------------------------------------
##
## accessors  [ should any be moved to svclasses? ]
##
##--------------------------------------------------
setMethod("chromosome", "GAlignments", function(object) as.character(seqnames(object)))
setMethod("first", "GAlignmentsList", function(x) x[[1]])
setMethod("last", "GAlignmentsList", function(x) x[[2]])

setReplaceMethod("first", "GAlignmentPairs", function(x, value){
  x@first <- value
  x
})
setReplaceMethod("last", "GAlignmentPairs", function(x, value){
  x@last <- value
  x
})
##
##--------------------------------------------------


##--------------------------------------------------
##
## Rearrangment functions
##
##--------------------------------------------------

bothTagsMapToCluster <- function(gpairs, tag_cluster){
  is_true <- overlapsAny(first(gpairs), tag_cluster) &
    overlapsAny(last(gpairs), tag_cluster)
}

annotateRPsWithLinkId <- function(gpairs, tag_cluster){
  F <- first(gpairs)
  L <- last(gpairs)
  ##
  ## Remove all read pairs for which the first tag and the last tag
  ## map to the same cluster
  ##
  hitsF <- findOverlaps(F, tag_cluster)
  hitsL <- findOverlaps(L, tag_cluster)
  if(is.null(names(tag_cluster))){
    nms <- seq_len(length(tag_cluster))
  } else nms <-  names(tag_cluster)
  mcols(F)$rpid <- nms[subjectHits(hitsF)]
  mcols(L)$rpid <- nms[subjectHits(hitsL)]
  first(gpairs) <- F
  last(gpairs) <- L
  gpairs
}

rpclustNamesMatrix <- function(nms){
  nms2 <- do.call(rbind, strsplit(nms[!duplicated(nms)], "-"))
  colnames(nms2) <- c("rpidF", "rpidL")
  nms2
}

linkClustersByReadPairs <- function(clustered_tags, gpairs){
  linked_cluster_ids <- rpclustNamesMatrix(names(gpairs))
  linked_clusters <- clustered_tags[linked_cluster_ids[, 1]]
  linked_clusters$linked.to <- clustered_tags[linked_cluster_ids[, 2]]
  names(linked_clusters) <- paste(linked_cluster_ids[, 1], linked_cluster_ids[, 2], sep="-")
  linked_clusters
}

atLeastOneTagMapsToCluster <- function(gpairs, tag_cluster){
  is_true <- overlapsAny(first(gpairs), tag_cluster) |
    overlapsAny(last(gpairs), tag_cluster)
}

gapToGRanges <- function(gpairs, linked.granges){
  atleast.one <- atLeastOneTagMapsToCluster(gpairs, linked.granges)
  ##checkIdentical(sum(atleast.one), 2949L)
  irp.nonlinking <- gpairs[ atleast.one ] ## why called "nonlinking"
  ## convert the GAlignmentPairs for these tags to a GRanges object
  mcF <- mcols(first(irp.nonlinking))
  F <- granges(first(irp.nonlinking))
  mcL <- mcols(last(irp.nonlinking))
  L <- granges(last(irp.nonlinking))
  tags <- c(F, L)
  linked.to <- c(L, F)
  tags$linked.to <- linked.to
  ## because we've concatentated F and L as a long vector, we can
  ## remove tags that do not map (we will have it in the other
  ## orientation)
  tags <- subsetByOverlaps(tags, linked.granges)
  tags
}

mapTagsToCluster <- function(tags, linked.clusters){
  hits1 <- findOverlaps(tags, linked.clusters, select="first")
  hits2 <- findOverlaps(tags, linked.clusters$linked.to, select="first")
  linked.ids <- strsplit(names(linked.clusters), "-")
  linked1id <- sapply(linked.ids, "[", 1)
  linked2id <- sapply(linked.ids, "[", 2)
  i1 <- linked1id[hits1]
  i2 <- linked2id[hits2]
  ##length(unique(names(index))
  mat <- cbind(i1, i2)
  rs <- rowSums(is.na(mat))
  ## for all the rows with no NAs, verify, the id is the same
  mat2 <- mat[rs==0, ]
  if(!identical(mat2[, 1], mat2[, 2])) stop(" the linked bin id to which the tags map are not the same ")
  if(any(rs==2)) stop("some of the tags did not map to any cluster")
  mapping <- rep(NA, length(tags))
  mapping[!is.na(i1)] <- i1[!is.na(i1)]
  mapping[is.na(i1)] <- i2[is.na(i1)]
  mapping
}

rpclustNames <- function(gpairs){
  rpidF <- rpid1 <- mcols(first(gpairs))$rpid
  rpidL <- rpid2 <- mcols(last(gpairs))$rpid
  rpindex1 <- as.integer(gsub("rpclust", "", rpid1))
  rpindex2 <- as.integer(gsub("rpclust", "", rpid2))
  rpidF[rpindex2 < rpindex1] <- rpid2[rpindex2 < rpindex1]
  rpidL[rpindex2 < rpindex1] <- rpid1[rpindex2 < rpindex1]
  paste0(rpidF, "-", rpidL)
}

Rearrangement <- function(linkedBins=GRanges(linked.to=GRanges(),
                                             partition=integer()),
                          improper=.GAlignmentPairs(),
                          partitioning=integer(),
                          unlinked_clusters=GRanges(),
                          linked1id=character(),
                          linked2id=character(),
                          tags=GRanges(linked.to=GRanges()),
                          tag_map_to_linked_bin=character()){
  if(length(unlinked_clusters) > 0 && length(improper) > 0 ){
    if(is.null(names(unlinked_clusters))){
      names(unlinked_clusters) <- seq_along(unlinked_clusters)
    }
    is_linking <- bothTagsMapToCluster(improper, unlinked_clusters)
    irp <- improper[ is_linking ]
    irp <- annotateRPsWithLinkId(irp, unlinked_clusters)
    names(irp) <- rpclustNames(irp)
    irp.nonlinking <- improper
    rpnames <- names(irp)
    linkedBins <- linkClustersByReadPairs(unlinked_clusters, irp)
    ## Many of the clusers are not linked to any other cluster.
    ## Because the unlinked_clusters are non-overlapping, these can be
    ## safely removed.
    is_true <- unlinked_clusters %in% linkedBins |
      unlinked_clusters %in% linkedBins$linked.to
    unfiltered_clusters <- unlinked_clusters[ is_true ]
    ##
    ## the first partitioning of the improper read pairs is by those
    ## that provide the link between the 'linked bins'
    ##
    partition <- setNames(seq_len(length(linkedBins)), names(linkedBins))
    partitioning <- partition[rpnames]
    ##improper <- irp
    linked.ids <- strsplit(names(linkedBins), "-")
    linked1id <- sapply(linked.ids, "[", 1)
    linked2id <- sapply(linked.ids, "[", 2)
    ##
    ## Tags that do not link clusters that were selected according to
    ## above filters, but that might be important for deciding whether
    ## rearrangement is an artifact
    ##
    ## convert galignmentpairs to granges
    ## Require at least one tag to map to a cluster
    tags <- gapToGRanges(irp.nonlinking, unfiltered_clusters)
    ##
    ## While the unlinked_clusters are non-overlapping, the linkedBins
    ## can have duplicates.  This means that a single tag will map to
    ## multiple linkedBins.  Subsetting a linkedBins that is
    ## duplicated should pull out the same set of single tags
    ##
    mapping <- mapTagsToCluster(tags, linkedBins)
    improper <- irp
  } else{
    mapping <- character()
  }
  lt <- linkedBins$linked.to
  seqlevels(lt, pruning.mode="coarse") <- seqlevels(lt)
  seqinfo(lt) <- seqinfo(linkedBins)
  linkedBins$linked.to <- lt
  if(length(linkedBins) == 0){
    names(linkedBins) <- character()
  }
  new("Rearrangement", linkedBins=linkedBins,
      link1id=linked1id, link2id=linked2id,
      partitioning=partitioning,
      improper=improper, tags=tags,
      tag_map_to_linked_bin=mapping)
}

##--------------------------------------------------
##
##  Tools to find rearrangements
##
##--------------------------------------------------

unpairThenReduce <- function(gpairs, params){
  g <- as(c(first(gpairs), last(gpairs)), "GRanges")
  strand(g) <- factor("*", levels=c("+","-","*"))
  rg <- reduce(g, min.gapwidth=minGapWidth(params))
  rg
}

filterTagClusters <- function(clustered_tags, gpairs, params){
  ## Remove clusters with low tag counts
  counts <- (countOverlaps(clustered_tags, first(gpairs)) +
             countOverlaps(clustered_tags, last(gpairs)))
  clustered_tags[ counts >= minNumberTagsPerCluster(params) ]
}

clusterTags <- function(gpairs, params){
  unpaired <- unpairThenReduce(gpairs, params)
  unpaired <- unpaired[ width(unpaired) >= minClusterSize(params) &
                        width(unpaired) <= maxClusterSize(params) ]
  unpaired <- unpaired[ !duplicated(unpaired) ]
  names(unpaired) <- seq_len(length(unpaired))
  unpaired <- filterTagClusters(unpaired, gpairs, params)
  unpaired
}

.trimInvalidReads <- function(x){
  chromosome <- function(x) as.character(seqnames(x))
  ends <- setNames(end(first(x)), chromosome(first(x)))
  ends.info <- seqlengths(x)[names(ends)]
  is_valid1 <- ends <= ends.info
  ends2 <- setNames(end(last(x)), chromosome(last(x)))
  ends.info2 <- seqlengths(x)[names(ends2)]
  is_valid2 <- ends2 <= ends.info2
  is_valid <- is_valid1 & is_valid2
  x[is_valid]
}

setMethod("numberLinkingRP", "Rearrangement", function(object){
  nsupportingReads <- table(names(improper(object)))
  nsupportingReads <- nsupportingReads[names(object)]
  as.integer(nsupportingReads)
})


seqJunctionsInferredByPairedTags <- function(aln.file, bins, param){
  gp <- readRDS(aln.file)
  ##gp <- updateObject(gp)
  seqlevelsStyle(gp) <- "UCSC"
  gp <- keepSeqlevels(gp, seqlevels(bins), pruning.mode="coarse")
  gp <- .trimInvalidReads(gp)
  gp2 <- filterPairedReads(gp, bins, param)
  ctags <- clusterTags(gp2, param)
  ##
  ## Note, the Rearrangement constructor does additional work
  ##
  L <- Rearrangement(improper=gp2, unlinked_clusters=ctags)
  L2 <- L[ numberLinkingRP(L) >= minNumberTagsPerCluster(param) ]
  L2
}

mapq <- function(galp){
  f <- first(galp)
  l <- last(galp)
  cbind(mcols(f)$mapq, mcols(l)$mapq)
}

#' Constructs a Rearrangement object of unlinked tag-clusters
#'
#' Identifies genomic intervals containing clusters of improper reads
#'  and constructs a Rearrangement object containing the unlinked
#'  clusters and the improper read pairs.
#'
#' @details
#' Read pairs are selected with parameters set in the
#' \code{RearrangementParams} object (denoted \code{params}) as follows:
#'
#' 1. the distance between the first and last read for a given pair
#' must be at least \code{rp_separation(params)}
#'
#' 2. Duplicate read pairs are dropped
#'
#' 3. For reads passing (1) and (2), the total number of reads aligned
#' to each bin in \code{bins} is counted. We exclude bins for which
#' the number of aligned reads is less than
#' \code{minNumberTagsPerCluster(params)}.  We use the remaining bins
#' to subset the read pairs object -- keeping only read pairs for
#' which the first or the last read overlaps a bin. Steps 1-3 are
#' performed by the function \code{filterPairedReads}.
#'
#' 4. Clusters of improper reads that pass (1), (2), and (3) are
#'   identified. In particular, we apply the function
#'   \code{unpairThenReduce} that first uncouples the read pairs and
#'   creates a single \code{GRanges} object.  The \code{GRanges}
#'   object is then reduced with argument \code{min.gapwidth} set to
#'   \code{minGapWidth(params)} (default is 1kb).  The interval for
#'   each cluster is given by the reduced ranges.  We keep only those
#'   intervals that have a width of at least
#'   \code{minClusterSize(params)} (default 115bp) and no larger than
#'   \code{maxClusterSize(params)} (default 5000bp).  In addition, we
#'   keep only those intervals for which the number of reads belonging
#'   to the interval is at least
#'   \code{minNumberTagsPerCluster(params)} (default 5). Step 4 is
#'   performed by the function \code{clusterTags}.
#'
#' 5.  Given a set of unlinked clusters and improper read pairs
#'   (filtered by steps 1-4), the Rearrangement constructor does the
#'   following:
#'
#'   i. annotates the red pairs with a unique id for cluster membership
#'
#'   ii. links the clusters (called linkedBins) (\code{linkClustersByReadPairs})
#'
#'   iii.  partitions the improper read pairs according to whether
#'   they link two clusters (a read pair can belong to multiple paired
#'   clusters)
#'
#'   iv.  maps each tag to a cluster
#'
#' REFACTORING: Rearrangement should do nothing and should be able to
#'   construct an empty Rearrangement object if no data is provided.
#'   Move the functions that do step 5 out of the constructor.
#'
#' @examples
#'   data(pdata, package="trellis")
#'   rp <- RearrangementParams()
#'   ##
#'   ## The file of improper read pairs is large, so this is slow
#'   ##
#'   r <- seqJunctionsInferredByPairedTags2(pdata,
#'                                          param=rp)
#'   r
#'   ## improper read pairs that link the clusters
#'   head(improper(r))
#'   ## The linked tag cluster intervals
#'   linkedBins(r)
#' @export
#'
#' @param preprocess a list as created by \code{preprocessData}
#'
#' @param param A \code{RearrangementParams} object
seqJunctionsInferredByPairedTags2 <- function(preprocess, param){
  bins <- preprocess$bins
  gp <- preprocess$read_pairs[["improper"]]
  m <- mapq(gp)
  gp <- gp[rowSums(m < 30) == 0]
  gp <- keepSeqlevels(gp, seqlevels(bins), pruning.mode="coarse")
  gp2 <- filterPairedReads(gp, bins, param)
  ctags <- clusterTags(gp2, param)
  ##
  ## Note, the Rearrangement constructor does additional work
  ##
  L <- Rearrangement(improper=gp2, unlinked_clusters=ctags)
  L2 <- L[ numberLinkingRP(L) >= minNumberTagsPerCluster(param) ]
  L2
}

computeFractionLinkingTags <- function(object){
  length(improper(object))*2/length(tags(object))
}

adjustLeftBinSize <- function(F, L){
  chrom.F <- as.integer(factor(chromosome(F), levels=seqlevels(F)))
  chrom.L <- as.integer(factor(chromosome(L), levels=seqlevels(L)))
  if(all(chrom.F != chrom.L)){
    is_less <- chrom.F < chrom.L
  } else {
    is_less <- start(F) < start(L)
  }
  ends <- starts <- rep(NA, length(F))
  if(any(is_less)){
    starts[ is_less ] <- start(F)[ is_less ]
    ends[ is_less ] <- end(F)[ is_less ]
    chrom.left <- (chromosome(F)[is_less])[1]
  }
  if(!all(is_less)){
    starts[ !is_less ] <- start(L)[ !is_less ]
    ends[ !is_less ] <- end(L)[ !is_less ]
    chrom.left <- (chromosome(L)[ !is_less ])[1]
  }
  left.bin <- GRanges(chrom.left, IRanges(min(starts), max(ends)))
  left.bin
}

adjustRightBinSize <- function(F, L){
  chrom.F <- as.integer(factor(chromosome(F), levels=seqlevels(F)))
  chrom.L <- as.integer(factor(chromosome(L), levels=seqlevels(L)))
  if(all(chrom.F != chrom.L)){
    is_greater <- chrom.F < chrom.L
  } else {
    is_greater <- start(F) < start(L)
  }
  ends <- starts <- rep(NA, length(F))
  if(any(is_greater)){
    starts[ is_greater ] <- start(L)[ is_greater ]
    ends[ is_greater ] <- end(L)[ is_greater ]
    chrom.right <- (chromosome(L)[is_greater])[1]
  }
  if(!all(is_greater)){
    starts[ !is_greater ] <- start(F)[ !is_greater ]
    ends[ !is_greater ] <- end(F)[ !is_greater ]
    chrom.right <- (chromosome(F)[ !is_greater ])[1]
  }
  right.bin <- GRanges(chrom.right, IRanges(min(starts), max(ends)))
  right.bin
}

adjustBinSizeByLinkedTagPairs <- function(object){
  linked_tags <- improper(object)
  si <- seqinfo(linked_tags)
  if(length(object) > 1) stop("object must be of length one (only one rearrangement)")
  F <- first(linked_tags)
  L <- last(linked_tags)
  left.bin <- adjustLeftBinSize(F, L)
  right.bin <- adjustRightBinSize(F, L)
  left.bin$linked.to <- right.bin
  linked.bin <- left.bin
  names(linked.bin) <- names(linkedBins(object))
  seqlevels(linked.bin, pruning.mode="coarse") <- seqlevels(si)
  seqinfo(linked.bin) <- si
  linked.bin
}

isTranslocation <- function(aln){
  chr1 <- seqnames(first(aln))
  chr2 <- seqnames(last(aln))
  as.logical(chr1!=chr2)
}

R1strand <- function(gpairs) strand(first(gpairs))
R2strand <- function(gpairs) strand(last(gpairs))
R1lessR2 <- function(gpairs){
  is.trans <- isTranslocation(gpairs)
  is.r1less <- as.logical(start(first(gpairs)) < start(last(gpairs)))
  is.r1less[is.trans] <- NA
  is.r1less
}
R1pos <- function(gpairs) as.logical(R1strand(gpairs)=="+")
R2neg <- function(gpairs) as.logical(R2strand(gpairs)=="-")

typeGAPairs <- function(gpairs){
  extdir <- system.file("extdata", package="trellis")
  rear_type <- read.csv(file.path(extdir, "rearrangement_types.csv"), nrows=20,
                        stringsAsFactors=FALSE)
  ##rear_type <- rear_type2
  rtypes <- matrix(NA, length(gpairs), nrow(rear_type))
  is.r1pos <- R1pos(gpairs)
  is.r1neg <- !is.r1pos
  is.r2neg <- R2neg(gpairs)
  is.r2pos <- !is.r2neg
  is.r1less <- R1lessR2(gpairs)
  is.r2less <- !is.r1less
  is.trans <- isTranslocation(gpairs)
  not.trans <- !is.trans
  is.aberrant.sep <- aberrantSep(gpairs)
  not.aberrant.sep <- !is.aberrant.sep
  ##colnames(rtypes) <- rear_type$type
  colnames(rtypes) <- rear_type$name
  rtypes[, 1] <- is.r1pos & is.r2neg & is.r1less & not.trans & not.aberrant.sep
  rtypes[, 2] <- is.r1neg & is.r2pos & is.r2less & not.trans & not.aberrant.sep
  rtypes[, 3] <- is.r1pos & is.r2neg & is.r1less & not.trans & is.aberrant.sep
  rtypes[, 4] <- is.r1neg & is.r2pos & is.r2less & not.trans & is.aberrant.sep
  rtypes[, 5] <- is.r1pos & is.r2neg & is.r1less & not.trans & not.aberrant.sep
  rtypes[, 6] <- is.r1pos & is.r2neg & is.r2less & not.trans
  rtypes[, 7] <- is.r1neg & is.r2pos & is.r2less & not.trans & not.aberrant.sep
  rtypes[, 8] <- is.r1neg & is.r2pos & is.r1less & not.trans
  rtypes[, 9] <-  is.r1pos & is.r2neg & is.trans
  rtypes[, 10] <- is.r1neg & is.r2pos & is.trans
  rtypes[, 11] <- is.r1pos & is.r2neg & is.trans
  rtypes[, 12] <- is.r1neg & is.r2pos & is.trans
  ## inversions
  rtypes[, 13] <- is.r1pos & is.r2pos & is.r1less & not.trans
  rtypes[, 14] <- is.r1pos & is.r2pos & is.r2less & not.trans
  rtypes[, 15] <- is.r1neg & is.r2neg & is.r1less & not.trans
  rtypes[, 16] <- is.r1neg & is.r2neg & is.r2less & not.trans
  rtypes[, 17] <- is.r1pos & is.r2pos &  is.trans
  rtypes[, 18] <- is.r1pos & is.r2pos &  is.trans
  rtypes[, 19] <- is.r1neg & is.r2neg &  is.trans
  rtypes[, 20] <- is.r1neg & is.r2neg &  is.trans
  col.indices <- split(seq_len(ncol(rtypes)), rear_type$compatible)
  y <- vector("list", length(col.indices))
  ##y <- foreach(j = col.indices, .combine="cbind") %do% {
  for(i in seq_along(col.indices)){
    j <- col.indices[[i]]
    y[[i]] <- rowSums(rtypes[, j, drop=FALSE])
  }
  y <- do.call(cbind, y)
  nms <- sapply(split(rear_type$name, rear_type$compatible), function(x) paste(x, collapse=","))
  colnames(y) <- nms
  y
}

#' Determine the type of rearrangement supported by each improper read pair
#'
#' @details
#' 
#' Rearrangements are typed by the following criteria for read 1 (R1)
#'   and read 2 (R2):
#' 
#' \enumerate{
#' 
#'   \item strand of R1 and R2
#'   \item orientation of R1 and R2 (e.g., R1 < R2)
#'   \item whether R1 and R2 are on different chromosomes
#'   \item whether R1 and R2 have an aberrant separation
#' }
#'
#' \strong{Deletions:} A deletion results in a single new sequence
#'   junction.  There are two orientations of a R1 and R2 that support
#'   this junction: R1+ < R2- and R1- > R2+.  The former R1+ < R2- is
#'   the sequenced + strand DNA fragment, while the latter R1- > R2+
#'   is the sequenced - strand.  These orientations are the expected
#'   orientations in unrearranged genomes, except that the distance
#'   between R1 and its mate is larger than expected (> 10kb).
#'
#' \strong{Amplicons:}
#' For purpose of downstream analyses, we only type intra-chromosomal
#'   amplicons that are replicated at the same site in the genome and
#'   without any changes to the strand orientation.  Illustration:
#'
#' \preformatted{
#'  Reference:
#'         1  2  3     4 5 6    7 8 9
#'   +5' ------------|--------|-------------
#'   -3' ------------|--------|-------------
#' 
#'  Tumor
#'         1  2  3     4 5 6    4 5 6   7 8 9
#'   +5' ------------|--------X-------|---------
#'   -3' ------------|--------X-------|---------
#'  }
#' Note that only 'X' is a new sequence junction not seen in the
#'   reference.  Characteristics of 'X' are R1+ > R2- (+ fragment) and
#'   R1- < R2+ (- fragment).  The amplicon could also insert further
#'   downstream:
#'
#' \preformatted{
#'  Tumor 
#'         1  2  3     4 5 6    7 8 9   4 5 6
#'   +5' ------------|--------|-------X---------
#'   -3' ------------|--------|-------X---------
#'
#'  or further upstream:
#' 
#'  Tumor 
#'         4  5  6     1 2 3    4 5 6   7 8 9
#'   +5' ------------X--------|-------|---------
#'   -3' ------------X--------|-------|---------
#'
#' }
#' 
#'   In each case, the rearrangement junctions given by 'X' would be
#'   identified by R1+ > R2 - and R1- < R2+.  
#'
#' \strong{Inter-chromosomal translocations:} For translocations,
#'   positional orientation is determined by chromosome number and not
#'   genomic position.  WLOG, we consider chr1 to be less than chr2
#'   and chr22 to be less that chrX.  In the case of an unbalanced
#'   translocation, a single rearrangement junction will be identified
#'   by the improper read pairs.
#'
#' \preformatted{
#' Reference   (-- chr1, == chr2)
#'   
#'      chr1  1 2 3   4 5 6   chr2  20 21 22   23  24  25
#' 5' + ------------|------   ===============|==============
#' 3' - ------------|------   ===============|==============
#' 
#' Tumor
#' 
#'       chr1  1 2 3   23 24 25 
#' 5' + ------------|========== 
#' 3' - ------------|==========
#' 
#' }
#'
#' Read orientations R1+ < R2- and R1- > R2+ are consistent with the
#'   fusion of the two positive strands and the two negative strands,
#'   respectively.  If the translocation is balanced, we would also
#'   observe
#'
#' \preformatted{
#' 
#' Tumor
#'      chr2  20 21 22   4 5 6
#' 5' + ==============|--------
#' 3' - ==============|--------
#'
#' }
#'
#' with R1+ > R2- and R1- < R2+.  Hence, there are four distinct read
#'   pair orientations for a balanced tranlocation. For the purpose of
#'   assessing gene fusions, we treat all translocations as if they
#'   are balanced even if we do not observe all 4 possible
#'   orientations.  An inversion translocation is typed and analyzed
#'   as an inversion.
#'
#'
#' \strong{Inversions:}
#'
#' An intrachromosomal inversion: 
#'
#' \preformatted{
#' 
#' Reference:  (-- positive strand,  == negative strand)
#' 
#'       1 2 3   4 5 6     7 8 9
#' 5'+  ------>|-------->|------->
#' 3'-  <======|<========|<=======
#'
#' Tumor (inversion)
#' 
#'       1 2 3   6 5 4    7 8 9
#' 5'+  ------>X=======>X------->
#' 3'-  <======X<-------X<=======
#'
#' }
#'
#' Again, X's denote the new sequence junctions formed as a result of
#'   the inversion.  The left-most X's are supported by R1+ < R2+ (the
#'   top strand in the diagram) and R1+ > R1+ (bottom strand).  The
#'   right-most X's are supported by R1- < R2- (top strand) and R1- >
#'   R2- (bottom strand).
#'
#' An inversion can also involve a translocation.
#'
#' @note Type \code{amp1,amp3} would never be identified because the
#'   junction does not involve an aberrant separation between read
#'   pairs.
#' 
#' @return a \code{data.frame} with colnames 'type' and 'percent'
#' @export
#' @param object a \code{Rearrangement} 
rearrangementType <- function(object){
  imp <- improper(object)
  rtypes <- typeGAPairs(imp)
  mns <- pmin(colMeans(rtypes), 1)
  ##nms <- .rearrangement_types()
  nms <- colnames(rtypes)
  data.frame(type=nms[which.max(mns)],
             percent=max(mns), stringsAsFactors=FALSE)
}

#' Finds candidate somatic rearrangements
#'
#' This function identifies clusters of improper reads that are linked
#' by the mate information in paired read sequencing platforms such as
#' Illumina HiSeq.
#'
#' @details
#' 
#' All reads from improper read pairs  where mates are separated by
#' at least 10kb and both reads in pair are mapped are read from the
#' \code{AlignmentViews} object.  A cluster of reads (all involved in
#' improper pairs) is defined as follows:
#' \enumerate{
#' 
#'   \item genomic intervals demarcating improper read clusters are
#'   gotten by applying reduce to a GRanges representation of all
#'   improper reads
#'
#'   \item genomic intervals must be at least 115bp and no larger than
#'   5000bp (default settings)
#'
#'   \item each cluster must contain at least 5 reads
#'
#' }
#' 
#' Non-overlapping clusters that are linked by multiple improper read
#' pairs are suggestive of a rearrangement.  Linked tag clusters are
#' identified by the function \code{seqJunctionsInferredByPairedTags}.
#' The genomic intervals defined by the linked tag clusters (also
#' referred to as linked bins) are represented as a \code{GRanges}
#' object with a variable called \code{linked.to} in \code{mcols}.
#' The \code{linked.to} column is also a \code{GRanges} object. The
#' \code{GRanges} object of the linked clusters, the improper read
#' pairs supporting the link, and the set of all tags that map to
#' either linked genomic interval are encapsulated in a
#' \code{Rearrangement} object.  Statistics calculated on each
#' \code{Rearrangement} object include the fraction of all reads link
#' the two clusters (\code{fractionLinkingTags}), the types of
#' rearrangements supported (\code{rearrangementType}), the modal
#' rearrangement, and the percent of read pairs supporting the modal
#' rearrangement.  The collection of all linked clusters for a given
#' sample is represented as a \code{RearrangementList}.
#'
#' @seealso See \code{\link{seqJunctionsInferredByPairedTags2}} for
#'   additional details regarding the clustering of tags from improper
#'   pairs and the identification of linked tag clusters. See
#'   \code{\link{rearrangementType}} for the type of rearrangement
#'   supported by each read pair.  See \code{\link{preprocessData}} for constructing a list of elemented obtained from preprocessing.
#'
#' @examples
#' ## Load list of preprocessed data (see preprocessData)
#' data(pdata, package="trellis")
#' ## Parameters for finding candidate rearrangements
#' rparam <- RearrangementParams(min_number_tags_per_cluster=5,
#'                               rp_separation=10e3)
#' ## List of candidate rearrangements
#' rlist <- findCandidates2(pdata, rparam)
#' rlist
#' @export
#' @param preprocess A list of preprocessing data as constructed by \code{preprocessData}
#' @param rp A \code{RearrangementParams} object
findCandidates2 <- function(preprocess, rp=RearrangementParams()){
  ##file <- improperPaths(align_view)
  candidates <- suppressWarnings(
    seqJunctionsInferredByPairedTags2(preprocess, rp)
  )
  candidateList <- RearrangementList(candidates)
  candidateList <- type_each(candidateList)
  colData(candidateList)$modal_rearrangement <- modalRearrangement(candidateList)
  colData(candidateList)$percent_rearrangement <- percentRearrangement(candidateList)
  candidateList
}

#' Summary stats for rearrangement object
#'
#' @param r a \code{Rearrangement} object
type_rearrangement <- function(r){
  linkedBins(r) <- adjustBinSizeByLinkedTagPairs(r)
  fractionLinkingTags(r) <- computeFractionLinkingTags(r)
  rtypes <- rearrangementType(r)
  modalRearrangement(r) <- rtypes[["type"]]
  percentRearrangement(r) <- rtypes[["percent"]]
  return(r)
}

type_each <- function(rlist){
  for(i in seq_along(rlist)){
    rlist[[i]] <- type_rearrangement(rlist[[i]])
  }
  rlist
}

overlappingTranscripts <- function(r, build, maxgap=5000){
  if(missing(build)){
    build <- genome(improper(r))[[1]]
    if(is.na(build)){
      stop("As the genome build is not specified in the seqinfo, the build argument can not be missing.")
    }
  }
  pkgname <- paste0("svfilters.", build)
  data(transcripts, package=pkgname, envir=environment())
  g <- uncouple(linkedBins(r))
  tx <- subsetByOverlaps(transcripts, g, maxgap=maxgap)
  if(!all(overlapsAny(g, transcripts, maxgap=maxgap))){
    ## one or both regions do not overlap a transcript
    no.overlap <- !overlapsAny(g, transcripts, maxgap=maxgap)
    noncoding <- g[no.overlap]
    noncoding$tx_id <- ""
    noncoding$tx_name <- ""
    noncoding$gene_name <- paste0("noncoding", seq_len(sum(no.overlap)))
    noncoding$cancer_connection <- FALSE
    noncoding$biol_sign <- FALSE
    tx <- c(tx, noncoding)
  }
  hits <- findOverlaps(g, tx, maxgap=maxgap)
  tx <- tx[subjectHits(hits)]
  txlist <- split(tx, queryHits(hits))
  txlist <- lapply(txlist, function(g){
    if(length(g) == 1) {
      rgr <- granges(g)
      rgr$gene_name <- g$gene_name
      return(rgr)
    }
    rgr <- reduce(g, ignore.strand=TRUE, min.gapwidth=maxgap)
    rgr$gene_name <- paste(unique(g$gene_name), collapse=",")
    rgr
  })
  tx <- unlist(GRangesList(txlist))
  hits <- findOverlaps(g, tx, maxgap=maxgap, select="first")
  g$gene_name <- tx$gene_name[hits]
  g
}


check_splitreads <- function(r){
  sr <- splitReads(r)
  names(sr) <- paste0("sr", seq_along(sr))
  sr.list <- split(sr, sr$qname)
  lb <- uncouple(linkedBins(r))
  is.overlap <- lapply(sr.list, function(g, bins){
    is.overlap <- overlapsAny(bins, g, maxgap=50)
    if(all(is.overlap)){
      result <- rep(TRUE, length(g))
      return(result)
    }
    rep(FALSE, length(g))
  }, bins=lb)
  is.overlap <- unlist(is.overlap)
  is.overlap
}

cstrand <- function(g) as.character(strand(g))

typeRead <- function(gap){
  strands <- paste0(cstrand(first(gap)), cstrand(last(gap)))
}

#' Create a data.frame of rearranged read pairs and split reads supporting a rearrangement
#'
#' This function is useful for converting a Rearrangement object to a data.frame for subsequent visualization by ggplot.
#'
#' @param r a \code{Rearrangement} object
#' @param build character string indicating genome build (only hg19 and hg18 currently supported)
#' @param maxgap the allowable distance between a split or paired end read and a transcript for assessing whether coding regions
#' @examples
#'   extdata <- system.file("extdata", package="svbams")
#'   rfile <- file.path(extdata, "CGOV11T_1.bam.rds")
#'   rlist <- readRDS(rfile)
#'   r <- rlist[[1]]
#'   r2 <- fiveTo3Prime(r, "hg19")
#'   rearDataFrame(r2[[1]], "hg19") 
#' @export
rearDataFrame <- function(r, build, maxgap=5000){
  lb <- linkedBins(r)
  if(!"gene_name" %in% colnames(mcols(lb))){
    bins <- overlappingTranscripts(r, build, maxgap=maxgap)
    bins$gene_name <- make.unique(bins$gene_name)
    names(bins) <- c(r@link1id, r@link2id)
    bins$reverse <- FALSE ## by default
    bins2 <- bins[1]
    bins2$linked.to <- bins[2]
    linkedBins(r) <- bins2
  }
  .rear_DataFrame(r)
}

list_region_levels <- function(x){
  levels1 <- levels(x)
  levels2 <- paste0(c("5'-", "3'-"), levels1)
  list(levels1, levels2)
}

label_orientation <- function(x){
  levels.list <- list_region_levels(x)
  levels1 <- levels.list[[1]]
  levels2 <- levels.list[[2]]
  index <- which(x %in% levels1[1])
  regions <- rep(NA, length(x))
  regions[index] <- gsub(levels1[1], levels2[1], x[index])
  index <- which(x %in% levels1[2])
  regions[index] <- gsub(levels1[2], levels2[2], x[index])
  factor(regions, levels2)
}

make_levels_unique <- function(x, z){
  if(!identical(levels(x), levels(z))){
    return(z)
  }
  ## x 5-'label1 3'-label2
  ## z 5-'label1 3'-label2
  ## remap z levels
  ##   5-'label2 3'-label1
  x.levels <- levels(x)
  z <- as.character(z)
  z2 <- z
  z.levels <- c(gsub("3'", "5'", x.levels[2]),
                gsub("5'", "3'", x.levels[1]))
  z[z2==x.levels[1]] <- z.levels[1]
  z[z2==x.levels[2]] <- z.levels[2]
  z <- factor(z, levels=z.levels)
  z
}

#' Creates a data.frame with both possible 5 to 3 prime orientations
#'
#' @param maxgap integer. See \code{findOverlaps}
#' @param build character string -- only 'hg19' or 'hg18' supported
#' @param rlist a length-two list of orientations. Each element is an object of class \code{Rearrangement}
#' @export
rearDataFrameList <- function(rlist, build, maxgap=5000){
  n_split_reads <- lengths(lapply(rlist, splitReads))
  if(!all(n_split_reads > 0)) stop("Each rearrangement must have at least one split read")
  df1 <- rearDataFrame(rlist[[1]], build, maxgap)
  df1$region <- label_orientation(df1$region)
  df2 <- rearDataFrame(rlist[[2]], build, maxgap)
  region2 <- label_orientation(df2$region)
  df2$region <- make_levels_unique(df1$region, region2)
  df <- rbind(df1, df2)
  levels <- c(levels(df1$region),
              levels(df2$region))
  df$region <- factor(df$region, levels=levels)
  df
}

.rear_DataFrame <- function(r){
  bin1 <- linkedBins(r)
  bin2 <- linkedBins(r)$linked.to
  mcols(bin1) <- mcols(bin1)[, 1:2]
  bins <- setNames(c(bin1, bin2),
                   c(r@link1id, r@link2id))
  gap <- improper(r)
  names(gap) <- seq_len(length(gap))
  s <- strands(r)
  igr <- ga2gr(gap, is.improper=TRUE, use.mcols=TRUE)
  sr <- splitReads(r)
  regions <- bins$gene_name[match(igr$rpid, names(bins))]
  hits <- findOverlaps(sr, bins, ignore.strand=TRUE, maxgap=100)
  sr.regions <- bins$gene_name[subjectHits(hits)]
  if(length(regions) != length(igr)){
    stop("Unexpected inequality in GRanges lengths")
  }
  igr$region <- regions
  df <- as(igr, "data.frame") %>%
    as.tibble()
  df2 <- df[, 1:5]
  tagid <- df$tagid
  read <- paste0(df$read, df$strand)
  sr$qname <- sapply(strsplit(sr$qname, "\\."), "[", 1)
  ##is.reverse.strand <- any(sr$reverse)
  red.sr <- reduce(sr)
  ix <- findOverlaps(red.sr, c(granges(bin1), granges(bin2))) %>%
    subjectHits
  red.sr <- red.sr[ix]
  hits <- findOverlaps(sr, red.sr)
  ##
  ## what if reverse is not defined yet?
  ##
  is.rev <- split(sr$reverse[queryHits(hits)],
                  subjectHits(hits)) %>%
    map_lgl(any) %>%
    set_names(c("5p", "3p"))
  if(!is.rev["5p"]){
    junction_5p <- end(red.sr)[1]
  } else{
    junction_5p <- start(red.sr)[1]
  }
  if(!is.rev["3p"]){
    junction_3p <- start(red.sr)[2]  
  } else {
    junction_3p <- end(red.sr)[2]
  }
  df2$junction_5p <- sr$junction_5p <- junction_5p
  df2$junction_3p <- sr$junction_3p <- junction_3p
  if(!"qname" %in% colnames(df2))
    df2$qname <- NA
  sr.df <- as(sr, "data.frame")
  ##sr.df2 <- sr.df[, c(1:5, 13, 14)]
  sr.df2 <- sr.df[, colnames(df2)]
  ##browser()
  df2$read_type <- df2$read <- read
  ##df2$read_type <- paste0(substr(read, 1, 2), rep(s, 2))
  df2$reverse <- igr$reverse
  sr.df2$read <- "splitread"
  sr.df2$read_type <- "splitread"
  sr.df2$reverse <- sr$reverse
  df2$tagid <- tagid
  df2$junction_5p <- junction_5p
  df2$junction_3p <- junction_3p
  sr.df2$tagid <- as.numeric(factor(sr.df$qname)) +
    max(as.numeric(tagid))
  df3 <- rbind(df2, sr.df2)
  regions2 <- c(df$region, sr.regions)
  df3$region <- factor(regions2, levels=bins$gene_name)
  df3$tagid <- as.numeric(df3$tagid)
  ##
  ## if transcript is on reverse strand, we need to multiply the coordinates by a negative number 
  ##
  ##
  if(FALSE){
    df3$start[df3$reverse] <- -1*df3$start[df3$reverse]
    df3$end[df3$reverse] <- -1*df3$end[df3$reverse]
  }
  df3
}

readPairColors <- function(){
  cols <- brewer.pal(12, "Paired")
  cols <- cols[-c(7,8)]
  cols <- cols[1:8]
  cols <- c(cols, "black")
  names(cols) <- c("R1-+", "R2-+",
                   "R1+-", "R2+-",
                   "R1--", "R2--",
                   "R1++", "R2++",
                   "splitread")
  cols
}

readColors <- function(){
  cols <- brewer.pal(12, "Paired")
  cols <- cols[-c(7,8)]
  cols <- cols[1:5]
  ##cols <- c(cols, "black")
  names(cols) <- c("R1-", "R2+",
                   "R1+", "R2-",
                   ##"R1-", "R2-",
                   ##"R1+", "R2+",
                   "splitread")
  cols
}

## this is the old version
ggRearrange2 <- function(df){
  colors <- readPairColors()[unique(df$read_type)]
  nms <- names(readPairColors())
  df$read_type <- factor(df$read_type, levels=nms)
  n.levels <- length(levels(df$region))
  nc <- 2
  if(n.levels==4){
    nr <- 2;
  } else nr <- 1;
  absoluteNum <- function(x, ...) abs(x)
  read_type <- tagid <- NULL
  ggplot(df, aes(ymin=tagid-0.2, ymax=tagid+0.2,
                 xmin=start/1e6, xmax=end/1e6,
                 color=read_type,
                 fill=read_type, group=tagid)) +
    geom_rect() +
    ylab("read pair index") +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    scale_x_continuous(breaks=pretty_breaks(5), labels=absoluteNum) +
    xlab("Mb") +
    theme(axis.text.x=element_text(size=7, angle=45, hjust=1),
          axis.text.y=element_blank()) +
    facet_wrap(~region, scales="free_x", nrow=nr, ncol=nc)
}

peelLegend <- function(gg){
  if(!is(gg, "gg")) stop("object must be instance of class 'gg'")
  gtab <- ggplotGrob(gg)
  i <- which(sapply(gtab$grobs, function(x) x$name) == "guide-box")
  legend <- gtab$grobs[[i]]
  fig <- gg + theme(legend.position="none")
  gtab <- ggplotGrob(fig)
  list(gtab=gtab, legend=legend)
}

axis_limit5p <- function(df, basepairs){
  if(nrow(df) == 0) return()
   if(!df$reverse[1]){
     xlim1 <- c(min(df$start),
                df$junction_5p[1])
     padding <- basepairs - abs(diff(xlim1))
     xlim1[1] <- xlim1[1] - padding
  } else {
    xlim1 <- c(max(df$start),
               df$junction_5p[1])
    padding <- basepairs - abs(diff(xlim1))
    xlim1[1] <- xlim1[1] + padding
  }
  xlim1
}

axis_limit3p <- function(df, basepairs){
  if(!df$reverse[1]){
    xlim2 <- c(min(df$junction_3p),
               max(df$end))
    padding <- basepairs - abs(diff(xlim2))
    xlim2[2] <- xlim2[2] + padding
  } else {
    xlim2 <- c(df$junction_3p[1],
               min(df$end))
    padding <- basepairs - abs(diff(xlim2))
    xlim2[2] <- xlim2[2] - padding
  }
  xlim2
}

axis_limits <- function(df, basepairs){
  region <- NULL
  df1 <- filter(df, region==levels(region)[1])
  df2 <- filter(df, region==levels(region)[2])
  xlim1 <- axis_limit5p(df1, basepairs)
  xlim2 <- axis_limit3p(df2, basepairs)
  limits <- list(xlim1=xlim1, xlim2=xlim2)
  names(limits) <- levels(df$region)
  limits
}

axis_labels5p <- function(df, xlim1, num.ticks){
  absoluteNum <- function(x, ...) prettyNum(abs(x), big.mark=",")
  brks1 <- pretty(xlim1, n=num.ticks)
  if(!df$reverse[1]){
    brks1[length(brks1)] <- df$junction_5p[1]
    ## avoid overcrowding of labels
    delta <- brks1[length(brks1)] - brks1[(length(brks1)-1)]
    if(delta < 20){
      brks1 <- brks1[ (length(brks1)-1) ]
    }
  } else {
    brks1[1] <- df$junction_5p[1]
    delta <- abs(brks1[2] - brks1[1])
    if(delta < 20){
      brks1 <- brks1[-2]
    }
  }
  labs1 <- paste(absoluteNum(brks1), "bp")
  list(breaks=brks1, labels=labs1)
}

axis_labels3p <- function(df, xlim2, num.ticks){
  absoluteNum <- function(x, ...) prettyNum(abs(x), big.mark=",")
  brks2 <- pretty(xlim2, n=num.ticks)
  if(!df$reverse[1]){
    brks2[1] <- df$junction_3p[1]
    ## avoid crowded labels
    delta <- brks2[2] - brks2[1]
    if(delta < 20){
      brks2 <- brks2[-2]
    }
  } else {
    brks2[length(brks2)] <- df$junction_3p[1]
    ## avoid crowded labels
    delta <- abs(brks2[length(brks2)] - brks2[length(brks2)-1])
    if(delta < 20){
      brks2 <- brks2[-(length(brks2)-1)]
    }
  }
  labs2 <- paste(absoluteNum(brks2), "bp")
  list(breaks=brks2, labels=labs2)
}

.ggRearrange <- function(df, ylabel="Read pair index",
                         basepairs=400, num.ticks=5){
  colors <- readColors()[unique(df$read_type)]
  colors["splitread"] <- "black"
  nms <- names(readColors())
  df$read_type <- factor(df$read_type, levels=nms)
  region <- read_type <- tagid <- NULL
  df1 <- filter(df, region==levels(region)[1])
  df2 <- filter(df, region==levels(region)[2])
  limits <- axis_limits(df, basepairs)
  gene1 <- levels(df$region)[1]
  gene2 <- levels(df$region)[2]
  xlim1 <- limits[[gene1]]
  xlim2 <- limits[[gene2]]
  labs1 <- axis_labels5p(df1, xlim1, num.ticks)
  labs2 <- axis_labels3p(df2, xlim2, num.ticks)
  a <- ggplot(df1, aes(ymin=tagid-0.2,
                       ymax=tagid+0.2,
                       xmin=start,
                       xmax=end,
                       color=read_type,
                       fill=read_type, group=tagid)) +
    geom_rect() +
    ylab(ylabel) +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    scale_x_continuous(breaks=labs1[["breaks"]],
                       labels=labs1[["labels"]]) +
    coord_cartesian(xlim=xlim1) +
    xlab("") +
    theme(axis.text.x=element_text(size=7, angle=45, hjust=1),
          axis.text.y=element_blank(),
          plot.title=element_text(size=5)) +
    guides(color=FALSE, fill=FALSE) +
    geom_vline(xintercept=df$junction_5p[1], linetype="dashed") +
    ggtitle(paste0(df1$region[1], " (", df1$seqnames[1], ")"))
  if(df1$reverse[1]){
    a <- a + scale_x_reverse(breaks=labs1[["breaks"]],
                             labels=labs1[["labels"]])
  }
  b <- ggplot(df2, aes(ymin=tagid-0.2,
                       ymax=tagid+0.2,
                       xmin=start,
                       xmax=end,
                       color=read_type,
                       fill=read_type, group=tagid)) +
    geom_rect() +
    ylab("read pair index") +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    scale_x_continuous(breaks=labs2[["breaks"]],
                       labels=labs2[["labels"]]) +
    coord_cartesian(xlim=xlim2) +
    xlab("") +
    theme(axis.text.x=element_text(size=7, angle=45, hjust=1),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position="bottom",
          legend.direction="horizontal",
          plot.title=element_text(size=5)) +
    guides(color=FALSE, fill=FALSE) +
    geom_vline(xintercept=df$junction_3p[1], linetype="dashed") +
    ylab("") +
    ggtitle(paste0(df2$region[1], " (", df2$seqnames[1], ")"))
  if(df2$reverse[1]){
    b <- b + scale_x_reverse(breaks=labs2[["breaks"]],
                             labels=labs2[["labels"]])
  }
  ##
  ## plot both panels just to get the legend
  d <- ggplot(df, aes(ymin=tagid-0.2,
                      ymax=tagid+0.2,
                      xmin=start,
                      xmax=end,
                      color=read_type,
                      fill=read_type, group=tagid)) +
    geom_rect() +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    theme(legend.position="bottom", legend.direction="horizontal") +
    guides(color=guide_legend(title=""), fill=guide_legend(title=""))
  legend.grob <- peelLegend(d)[[2]]
  agrob <- ggplotGrob(a)
  bgrob <- ggplotGrob(b)
  ##legend.grob <- gg.objs[[2]]
  bgrob$widths <- agrob$widths
  list(`5p`=agrob,
       `3p`=bgrob,
       legend=legend.grob)
}


## plot sequences of split reads
.ggRearrange_sequences <- function(df, ylabel="Read pair index",
                                   basepairs=400, num.ticks=5){
  colors <- readColors()[unique(df$read_type)]
  colors["splitread"] <- "black"
  nms <- names(readColors())
  df$read_type <- factor(df$read_type, levels=nms)
  region <- read_type <- tagid <- NULL
  df1 <- filter(df, region==levels(region)[1])
  df2 <- filter(df, region==levels(region)[2])
  limits <- axis_limits(df, basepairs/2)
  gene1 <- levels(df$region)[1]
  gene2 <- levels(df$region)[2]
  xlim1 <- limits[[gene1]]
  xlim2 <- limits[[gene2]]
  labs1 <- axis_labels5p(df1, xlim1, num.ticks)
  labs2 <- axis_labels3p(df2, xlim2, num.ticks)
  seq2 <- seq1 <- NULL

  df1.seq <- filter(df1, !is.na(seq1))
  df2.seq <- filter(df2, !is.na(seq2))
  ##
  ## sequences are of variable lengths and have different starts and ends
  ## - need a data.frame of the form
  ## position qname dna_base
  position.list1 <- list()
  position.list2 <- list()
  base.list2 <- base.list1 <- list()
  for(i in seq_len(nrow(df1.seq))){
    position.list1[[i]] <- (df1.seq$start[i]+1):(df1.seq$end[i])
    position.list2[[i]] <- (df2.seq$start[i]+1):(df2.seq$end[i])
    tmp1 <- DNAString(df1.seq$seq1[[i]])
    tmp2 <- DNAString(df2.seq$seq2[[i]])
    base.list1[[i]] <- sapply(as.list(tmp1), as.character)
    base.list2[[i]] <- sapply(as.list(tmp2), as.character)
  }
  ##
  ##  match sequences to the dna coordinates based on the lenths
  ##
  sizes1 <- lengths(position.list1) 
  sizes2 <- lengths(base.list1)
  sizes3 <- lengths(base.list2)
  if(all(sizes1 == sizes2)){
    tab1 <- tibble(position=unlist(position.list1),
                   dna_base=unlist(base.list1))
    tab1$qname <- rep(df1.seq$qname, lengths(position.list1))
    tab1$tagid <- rep(df1.seq$tagid, lengths(position.list1))
    tab2 <- tibble(position=unlist(position.list2),
                   dna_base=unlist(base.list2))
    tab2$qname <- rep(df2.seq$qname, lengths(position.list2))
    tab2$tagid <- rep(df2.seq$tagid, lengths(position.list2))
    tab <- bind_rows(tab1, tab2)
  }
  if(all(sizes1 == sizes3)){
    tab1 <- tibble(position=unlist(position.list1),
                   dna_base=unlist(base.list2))
    tab1$qname <- rep(df1.seq$qname, lengths(position.list1))
     tab1$tagid <- rep(df1.seq$tagid, lengths(position.list1))
    tab2 <- tibble(position=unlist(position.list2),
                   dna_base=unlist(base.list1))
    tab2$qname <- rep(df2.seq$qname, lengths(position.list2))
    tab2$tagid <- rep(df2.seq$tagid, lengths(position.list2))
    tab <- bind_rows(tab1, tab2)
  }
  dna_base <- position <- NULL
  a <- ggplot(df1, aes(ymin=tagid-0.2,
                       ymax=tagid+0.2,
                       xmin=start,
                       xmax=end,
                       color=read_type,
                       fill=read_type, group=tagid)) +
    geom_rect() +
    geom_text(data=tab1,
              aes(x=position, y=tagid, label=dna_base),
              color="white",
              inherit.aes=FALSE,
              size=1) +
    ylab(ylabel) +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    scale_x_continuous(breaks=labs1[["breaks"]],
                       labels=labs1[["labels"]]) +
    coord_cartesian(xlim=xlim1) +
    xlab("") +
    theme(axis.text.x=element_text(size=7, angle=45, hjust=1),
          axis.text.y=element_blank(),
          plot.title=element_text(size=5)) +
    guides(color=FALSE, fill=FALSE) +
    geom_vline(xintercept=df$junction_5p[1], linetype="dashed") +
    ggtitle(paste0(df1$region[1], " (", df1$seqnames[1], ")"))
  if(df1$reverse[1]){
    a <- a + scale_x_reverse(breaks=labs1[["breaks"]],
                             labels=labs1[["labels"]])
  }
  b <- ggplot(df2, aes(ymin=tagid-0.2,
                       ymax=tagid+0.2,
                       xmin=start,
                       xmax=end,
                       color=read_type,
                       fill=read_type, group=tagid)) +
    geom_rect() +
    ylab("read pair index") +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    scale_x_continuous(breaks=labs2[["breaks"]],
                       labels=labs2[["labels"]]) +
    coord_cartesian(xlim=xlim2) +
    xlab("") +
    theme(axis.text.x=element_text(size=7, angle=45, hjust=1),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position="bottom",
          legend.direction="horizontal",
          plot.title=element_text(size=5)) +
    guides(color=FALSE, fill=FALSE) +
    geom_vline(xintercept=df$junction_3p[1], linetype="dashed") +
    ylab("") +
    ggtitle(paste0(df2$region[1], " (", df2$seqnames[1], ")"))
  if(df2$reverse[1]){
    b <- b + scale_x_reverse(breaks=labs2[["breaks"]],
                             labels=labs2[["labels"]])
  }
  ##
  ## plot both panels just to get the legend
  d <- ggplot(df, aes(ymin=tagid-0.2,
                      ymax=tagid+0.2,
                      xmin=start,
                      xmax=end,
                      color=read_type,
                      fill=read_type, group=tagid)) +
    geom_rect() +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    theme(legend.position="bottom", legend.direction="horizontal") +
    guides(color=guide_legend(title=""), fill=guide_legend(title=""))
  legend.grob <- peelLegend(d)[[2]]
  agrob <- ggplotGrob(a)
  bgrob <- ggplotGrob(b)
  ##legend.grob <- gg.objs[[2]]
  bgrob$widths <- agrob$widths
  list(`5p`=agrob,
       `3p`=bgrob,
       legend=legend.grob)
}



ggRearrangeLegend <- function(){
  read_type <- NULL
  cols <- readColors()
  cols[["splitread"]] <- "black"
  df <- data.frame(start=seq_along(cols),
                   end=seq_along(cols) + 1,
                   read_type=names(cols))
  fig <- ggplot(df, aes(xmin=start, xmax=end,
                        ymin=start, ymax=end,
                        color=read_type,
                        fill=read_type)) +
    geom_rect() +
    scale_fill_manual(values=cols) +
    scale_color_manual(values=cols) +
    guides(color=guide_legend(title=""),
           fill=guide_legend(title="")) +
    theme(legend.direction="horizontal")
  gobj <- peelLegend(fig)
  gobj[[2]]
}

#' ggplot wrapper for plotting reads supporting a rearrangement
#'
#' @param df a \code{data.frame} as created by \code{rearDataFrame}
#' @param ylab y-axis label
#' @param basepairs integer specifying size of window in basepairs
#' @param num.ticks integer specifying number of x-axis ticks
#' @seealso \code{\link{rearDataFrame}}
#' @examples
#'   data(rlist)
#'   extdata <- system.file("extdata", package="svbams")
#'   unmap.file <- file.path(extdata, "blat_unmapped.txt")
#'   blat_unmap <- readBlat(unmap.file)
#'   split_reads <- rearrangedReads(linkedBins(rlist),
#'                                  blat_unmap, 500)
#'   splitReads(rlist) <- split_reads
#'   rlist2 <- fiveTo3List(rlist, build="hg19")
#'   r <- rlist2[[1]]
#'   df <- rearDataFrame(r, "hg19")
#'   \dontrun{
#'     ggRearrange(df)
#'   }
#' @export
ggRearrange <- function(df, ylab="Read pair index",
                         basepairs=400, num.ticks=5){
  . <- NULL
  grobs <- .ggRearrange(df, ylabel=ylab,  basepairs, num.ticks)
  widths <- c(0.5, 0.5) %>%
    "/"(sum(.)) %>%
    unit(., "npc")
  heights <- c(0.95, 0.05) %>%
    "/"(sum(.)) %>%
    unit(., "npc")
  mat <- matrix(c(1, 2,
                  3, 3), byrow=TRUE, ncol=2, nrow=2)
  agrob <- grobs[["5p"]]
  bgrob <- grobs[["3p"]]
  legend.grob <- grobs[["legend"]]
  gobj <- grid.arrange(agrob, bgrob,
                       legend.grob,
                       layout_matrix=mat,
                       widths=widths,
                       heights=heights)
  list(arranged.grobs=gobj,
       grobs=grobs)
}


ggRearrangeSequences <- function(df, ylab="Read pair index",
                                 basepairs=400, num.ticks=5){
  . <- NULL
  grobs <- .ggRearrange_sequences(df, ylabel=ylab,  basepairs, num.ticks)
  widths <- c(0.5, 0.5) %>%
    "/"(sum(.)) %>%
    unit(., "npc")
  heights <- c(0.95, 0.05) %>%
    "/"(sum(.)) %>%
    unit(., "npc")
  mat <- matrix(c(1, 2,
                  3, 3), byrow=TRUE, ncol=2, nrow=2)
  agrob <- grobs[["5p"]]
  bgrob <- grobs[["3p"]]
  legend.grob <- grobs[["legend"]]
  gobj <- grid.arrange(agrob, bgrob,
                       legend.grob,
                       layout_matrix=mat,
                       widths=widths,
                       heights=heights)
  list(arranged.grobs=gobj,
       grobs=grobs)
}

strands <- function(r){
  gap <- improper(r)
  paste0(cstrand(first(gap)), cstrand(last(gap)))
}

isInversion <- function(r){
  s <- strands(r)
  any(s %in% c("++", "--"))
}

split_reads_order <- function(r, build, maxgap){
  bins <- overlappingTranscripts(r, build, maxgap=maxgap)
  bins$gene_name <- make.unique(bins$gene_name)
  names(bins) <- c(r@link1id, r@link2id)
  ##bins <- uncouple(linkedBins(r))
  names(bins) <- c(r@link1id, r@link2id)
  ## easiest might be to organize by split read
  sr <- splitReads(r)
  names(sr) <- paste0("sr", seq_along(sr))
  bins2 <- bins[subjectHits(findOverlaps(sr, bins, maxgap=50))]
  s <- strands(r)
  d1 <- abs(end(bins2) - end(sr))
  d2 <- abs(start(sr) - start(bins2))
  ratio <- d1/d2
  if(length(ratio) > 2){
    km <- kmeans(ratio, centers=2)$cluster
    means <- sapply(split(ratio, km), mean)
    ix <- which.min(means)
    is.fiveprime <- km[seq_along(sr)] == ix
  } else {
    ix <- which.min(ratio)
    is.fiveprime <- seq_along(sr) == ix
  }
  sr.fiveprime <- sr[is.fiveprime]
  no.regions <- length(reduce(sr.fiveprime))
  if(no.regions > 1){
    stop("Didn't expect to reach here 1")
  }
  five.prime.region <- overlapsAny(bins, sr[is.fiveprime])
  if(identical(five.prime.region, c(TRUE, FALSE))){
    ## linked bins already ordered 5' to 3'
    bins$reverse <- FALSE
    bins2 <- bins[1]
    bins2$linked.to <- bins[2]
    linkedBins(r) <- bins2
    return(r)
  }
  if(!identical(five.prime.region, c(FALSE, TRUE))){
    stop("Didn't expect to reach here 2")
  }
  r2 <- r
  bins$reverse <- FALSE
  bins2 <- bins[2]
  bins2$linked.to <- bins[1]
  ##lb <- linkedTo(r)
  ##lb$linked.to <- granges(linkedBins(r))
  link1id <- r@link2id
  link2id <- r@link1id
  names(bins2) <- paste(link1id, link2id, sep="-")
  linkedBins(r2) <- bins2
  r2@link1id <- link1id
  r2@link2id <- link2id
  splitReads(r2)$reverse <- FALSE
  ga <- improper(r2)
  mcols(first(ga))$reverse <- FALSE
  mcols(last(ga))$reverse <- FALSE
  r2@improper <- ga
  r2
}

list_orientations <- function(r1, r2){
  orientations <- list(r1, r2)
  names(orientations) <- c(names(linkedBins(r1)),
                           names(linkedBins(r2)))
  orientations
}



.type_rear <- function(object){
  ss <- strands(object)
  ss[ ss=="-+" ] <- "+-"
  ss <- unique(ss)
  ss <- paste(ss, collapse=",")
  ss
}

#' @aliases type,Rearrangement-method
#' @rdname Rearrangement-class
setMethod("type", "Rearrangement", function(object){
  .type_rear(object)
})

#' @rdname Rearrangement-class
#' @aliases type,RearrangementList-method
setMethod("type", "RearrangementList", function(object){
  x <- sapply(object, type)
  x
})

#' Ad-hoc assessment of the complexity of a rearrangement using the strand of the supporting read pairs
#'
#' For most rearrangements, the strand of the 5-prime and 3-prime reads belonging to a rearranged read pair will be consistent. For example, all rearranged read pairs are +/-, indicating a fusion on the positive strand (R1+, R2-) or negative strand (R1-, R2+).
#'
#' @param x a \code{Rearrangement} or \code{RearrangementList}
## @examples
##   extdata <- system.file("extdata", package="svbams")
##   r <- readRDS(file.path(extdata, "cgov1t_complex_rearrangement.rds"))
##   s <- type(r)
##   print(s)
##   isComplex(r)
isComplex <- function(x){
  ss <- type(x)
  xx <- strsplit(ss, ",")
  elementNROWS(xx) > 1
}


#' Put the linked genomic intervals in 5-prime to 3-prime order with respect to the rearranged genome
#'
#'
#' @param r a \code{Rearrangement} object
#' @param build string providing build of reference genome ('hg18', 'hg19', or 'hg38')
#' @param maxgap maximum distance between read and transcript for transcript to be considered overlapping
#' @seealso \code{\link{rearDataFrame}} \code{\link{ggRearrange}}
#' @examples
#'   extdata <- system.file("extdata", package="svbams")
#'   rfile <- file.path(extdata, "CGOV11T_1.bam.rds")
#'   rlist <- readRDS(rfile)
#'   r <- rlist[[1]]
#'   r2 <- fiveTo3Prime(r, "hg19")
#'   df <- rearDataFrame(r2[[1]], "hg19")
#'   ggRearrange(df)
#' @export
fiveTo3Prime <- function(r, build, maxgap=5000){
  s <- strands(r)
  tab <- sort(table(s), decreasing=TRUE)
  ## if both ++ and -- strands are observed
  ##  only the most frequent strands will be evaluated
  s <- names(tab)[1]
  orientations <- NULL
  if(all(s=="+-" | s=="-+")){
    ##r1 <- split_reads_order(r, build, maxgap)
    r1 <- posNeg(r, build, maxgap)
    r2 <- negPos(r, build, maxgap)
    orientations <- list_orientations(r1, r2)
  }
  if(all(s=="--")){
    r1 <- negativeInversion1(r, build, maxgap)
    r2 <- negativeInversion2(r, build, maxgap)
    orientations <- list_orientations(r1, r2)
  }
  if(all(s=="++")){
    r1 <- positiveInversion1(r, build, maxgap)
    r2 <- positiveInversion2(r, build, maxgap)
    orientations <- list_orientations(r1, r2)
  }
  return(orientations)
}

## This function is used when R1 and its mate align to the opposite strand
##
## posNeg refers to R1 on positive strand and R1's mate on the negative strand, or R1 on negative strand and R2 on positive strand
## - assumes R1 on positive is in bin 1.
##
posNeg <- function(r, build, maxgap=5000){
  bins <- overlappingTranscripts(r, build, maxgap=maxgap)
  bins$gene_name <- make.unique(bins$gene_name)
  names(bins) <- c(r@link1id, r@link2id)
  ##
  ## the improper pairs have an rpid column that indicates to which bin the read is linked
  ga <- improper(r)
  s <- strands(r)
  r1 <- first(ga)
  r2 <- last(ga)
  r1.rpid <- mcols(r1)$rpid
  r2.rpid <- mcols(r2)$rpid
  sr <- splitReads(r)
  if(length(sr) == 0) {
    message("no split reads")
    return(r)
  }
  sr <- subsetByOverlaps(sr, bins, maxgap=50)
  if (length(sr) == 0) {
    message(paste0("No split reads are within 50bp of improper read pairs for rearrangement: ", names(linkedBins(r))))
    return(r)
  }
  hits <- findOverlaps(sr, bins, maxgap=50)
  hits <- hits[!duplicated(queryHits(hits))]
  j <- subjectHits(hits)
  i <- queryHits(hits)
  sr$rpid <- NA
  sr$rpid[i] <- names(bins)[j]
  ##r1pos.is.bin1 <- r1.rpid == names(bins)[1] & cstrand(r1) == "+"
  ## if r1 is positive and in second bin, reverse orientation for both r1 and r2
  is.rev <- any(cstrand(r1)=="+" & r1.rpid == names(bins)[2])
  mcols(r2)$reverse <- mcols(r1)$reverse <- is.rev
  ## reverse all or reverse none
  sr$reverse <- is.rev[1]
  first(ga) <- r1
  last(ga) <- r2
  r@improper <- ga
  bins$reverse <- is.rev[1]
  bins$reverse[1] <- is.rev[1]
  bins2 <- bins[1]
  bins2$linked.to <- bins[2]
  names(bins2) <- paste(r@link1id, r@link2id, sep="-")
  linkedBins(r) <- bins2
  splitReads(r) <- sr
  r
}

negPos <- function(r, build, maxgap=5000){
  r <- swap_bin_order(r)
  r <- posNeg(r, build, maxgap)
  r
}

negativeInversion1 <- function(r, build, maxgap=5000){
  ##     multiple orientations
  ## *  i.  r1-- -> r5' and r2-- -> 3'  *
  ##   ii.  r1-- -> 3'  and r2-- -> r5'
  bins <- overlappingTranscripts(r, build, maxgap=maxgap)
  bins$gene_name <- make.unique(bins$gene_name)
  names(bins) <- c(r@link1id, r@link2id)
  ##
  ## the improper pairs have an rpid column that indicates to which bin the read is linked
  ga <- improper(r)
  s <- strands(r)
  r1 <- first(ga)
  r2 <- last(ga)
  r1.rpid <- mcols(r1)$rpid
  r2.rpid <- mcols(r2)$rpid
  sr <- splitReads(r)
  sr <- subsetByOverlaps(sr, bins, maxgap=50)
  sr$rpid <- names(bins)[subjectHits(findOverlaps(sr, bins, maxgap=50))]
  ##
  ##
  ## approx. half of the r1-- will be bin[1], the other have bin[2]
  ## here, we do i.
  r1.is.bin1 <- r1.rpid == names(bins)[1]
  mcols(r2)$reverse <- mcols(r1)$reverse <- FALSE
  mcols(r1)$reverse[r1.is.bin1] <- TRUE
  mcols(r2)$reverse[!r1.is.bin1] <- TRUE
  sr$reverse <- FALSE
  sr$reverse[sr$rpid==names(bins)[1]] <- TRUE
  first(ga) <- r1
  last(ga) <- r2
  r@improper <- ga
  bins$reverse <- FALSE
  bins$reverse[1] <- TRUE
  bins2 <- bins[1]
  bins2$linked.to <- bins[2]
  names(bins2) <- paste(r@link1id, r@link2id, sep="-")
  linkedBins(r) <- bins2
  splitReads(r) <- sr
  r
}

negativeInversion2 <- function(r, build, maxgap=5000){
  ##     multiple orientations
  ##   ii.  r1-- -> 3'  and r2-- -> r5'
  r <- swap_bin_order(r)
  r <- negativeInversion1(r, build, maxgap)
  r
}

positiveInversion1 <- function(r, build, maxgap=5000){
  ##     multiple orientations
  ## r1 -> 5', r2 -> reverse3'
  bins <- overlappingTranscripts(r, build, maxgap=maxgap)
  bins$gene_name <- make.unique(bins$gene_name)
  names(bins) <- c(r@link1id, r@link2id)
  ##
  ## the improper pairs have an rpid column that indicates to which bin the read is linked
  ga <- improper(r)
  s <- strands(r)
  r1 <- first(ga)
  r2 <- last(ga)
  r1.rpid <- mcols(r1)$rpid
  r2.rpid <- mcols(r2)$rpid
  sr <- splitReads(r)
  sr <- subsetByOverlaps(sr, bins, maxgap=50)
  sr$rpid <- names(bins)[subjectHits(findOverlaps(sr, bins, maxgap=50))]
  ##
  ##
  ## approx. half of the r1-- will be bin[1], the other have bin[2]
  ## here, we do i.
  r2.is.bin2 <- r2.rpid == names(bins)[2]
  mcols(r2)$reverse <- mcols(r1)$reverse <- FALSE
  mcols(r1)$reverse[!r2.is.bin2] <- TRUE
  mcols(r2)$reverse[r2.is.bin2] <- TRUE
  sr$reverse <- FALSE
  sr$reverse[sr$rpid==names(bins)[2]] <- TRUE
  first(ga) <- r1
  last(ga) <- r2
  r@improper <- ga
  bins$reverse <- FALSE
  bins$reverse[1] <- TRUE
  bins2 <- bins[1]
  bins2$linked.to <- bins[2]
  names(bins2) <- paste(r@link1id, r@link2id, sep="-")
  linkedBins(r) <- bins2
  splitReads(r) <- sr
  r
}

swap_bin_order <- function(r){
  bins <- linkedTo(r)
  bins2 <- linkedBins(r)
  drop <- colnames(mcols(bins2)) == "linked.to"
  mcols(bins2) <- mcols(bins2)[, !drop, drop=FALSE]
  bins$linked.to <- bins2
  l1id <- r@link2id
  l2id <- r@link1id
  r@link1id <- l1id
  r@link2id <- l2id
  names(bins) <- paste(l1id, l2id, sep="-")
  linkedBins(r) <- bins
  r
}

positiveInversion2 <- function(r, build, maxgap=5000){
  r <- swap_bin_order(r)
  r <- positiveInversion1(r, build, maxgap)
}

geneNames <- function(r){
  c(linkedBins(r)$gene_name, linkedTo(r)$gene_name)
}
