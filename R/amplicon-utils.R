#' @include help.R
NULL

setGeneric("zeroEndAnchors", function(object) standardGeneric("zeroEndAnchors"))
setGeneric("singleEndAnchors", function(object) standardGeneric("singleEndAnchors"))
setGeneric("bothEndAnchors", function(object) standardGeneric("bothEndAnchors"))

setClass("AnchoredReadPairs",
         representation(zero_end="GAlignmentPairs",
                        single_end="GAlignmentPairs",
                        both_ends="GAlignmentPairs"))

AnchoredReadPairs <- function(zero_end, single_end, both_ends=GRanges()){
  if(missing(zero_end)) new("AnchoredReadPairs")
  new("AnchoredReadPairs", zero_end=zero_end,
      single_end=single_end, both_ends=both_ends)
}

setMethod("zeroEndAnchors", "AnchoredReadPairs",
          function(object) object@zero_end)
setMethod("singleEndAnchors", "AnchoredReadPairs",
          function(object) object@single_end)
setMethod("bothEndAnchors", "AnchoredReadPairs",
          function(object) object@both_ends)


FilterEdgeParam <- function(freq=5, minimum_maxdist=50, bad_bins=GRanges()){
  list(freq=freq,
       minimum_maxdist=minimum_maxdist,
       bad_bins=bad_bins)
}

gaps0 <- function(x){
  g <- gaps(x)
  g <- g[strand(g)=="*"]
  g
}

hgnc <- function(g) g$hgnc

isDuplication <- function(object, minimum_foldchange=1){
  ## - must not already be called an amplicon
  ## - must be somatic
  ## - must have elevated copy number
  !isAmplicon(object) & (object$seg.mean >= minimum_foldchange) &
    !is.na(object$seg.mean) & isSomatic(object) & width(object) < 500e3
}


#' Remove the genomic intervals in query that overlap with intervals
#' in subject
#'
#'
#' @param query a \code{\linkS4class{BamViews}} or \code{\linkS4class{GRanges}} instance
#' @param subject a \code{\linkS4class{GRanges}} object
#' @param type see \code{\link{findOverlaps}}
#' @param ... Additional arguments passed to \code{findOverlaps}
#' @return Returns the \code{query} with intervals that overlap \code{subject} removed
#' @export
#' @seealso \code{\link[GenomicRanges]{findOverlaps}}
#' @rdname filterBy-methods
setGeneric("filterBy", function(query, subject, type="any", ...) standardGeneric("filterBy"))

#' @aliases filterBy,GRanges,GRanges-method
#' @rdname filterBy-methods
setMethod("filterBy", c("GRanges", "GRanges"),
          function(query, subject, type="any", ...){
            dropbins <- unique(queryHits(findOverlaps(query, subject, type=type, ...)))
            if(length(dropbins) > 0)
              query <- query[-dropbins]
            return(query)
          })

node2 <- function(name, sep="-") sapply(name, function(x) strsplit(x, sep)[[1]][2])
node1 <- function(name, sep="-") sapply(name, function(x) strsplit(x, sep)[[1]][1])

#' Lists germline filters and parameters for filtering amplicons
#'
#' 
#' @export
#'
#' @examples
#' library(svcnvs)
#' params <- ampliconParams("hg19")
#'
#' @return Returns a list of germline filters and some of the hard
#'   thresholds used in the amplicon analysis.
#'
#' @param build character string providing UCSC genome build
#'   (currently, must be hg19)
#'
#' @param AMP_THR numeric threshold for high copy amplicon
#'
#' @param LOW_THR numeric a lower threshold used when considering
#'   amplicons that are bridged by improper read pairs to a high copy
#'   amplicon
#'
#' @param border_size used to construct a query (GRanges object) for
#'   additional amplicons neighboring a focal amplicon. TODO: more
#'   detail needed.
#'
#' @param overhang An integer indicating how much to expand the
#'   germline filter on each size.  Using \code{overlapsAny}, a
#'   determination is made whether an amplicon is part-germline (any
#'   overlap) or no-germline (no overlap).  If the amplicon is
#'   completely within the extended germline filter, the amplicon is
#'   considered fully-germline.
#'
#' @param MIN_WIDTH length-one integer vector indicating the minimum size of amplicon (see also \code{joinNearGRanges})
#'
#' @param MIN_SEGMEAN_DIFF length-one numeric vector.  Adjacent segments whose means differ by less than this value are candidates for merging by \code{{joinNearGRanges}}
#'
#' @param min.gapwidth length-one numeric vector passed to the \code{reduce} method that merges adacent segments with comparable segment means
#'
#' @param maxgap length-one numeric vector. Passed to \code{findOverlaps} when
#'   evaluating whether the amplicon ranges overlap with other amplicons (see
#'   \code{linkNearAmplicons}}).
#'
#' @param frequency
#' @param minimum_maxdist
#' @param bad_bins a \code{GRanges} object of problematic bins
ampliconParams <- function(AMP_THR=log2(2.75),
                           LOW_THR=log2(1.75),
                           border_size=10e3,
                           overhang=25e3,
                           MIN_WIDTH=2000,
                           MIN_SEGMEAN_DIFF=0.05,
                           min.gapwidth=3000,
                           maxgap=500e3,
                           frequency=5,
                           minimum_maxdist=50,
                           minimum_count=5,
                           bad_bins=GRanges()){
  ##filters <- listGenomeFilters()
  edge <- FilterEdgeParam(freq=5,
                          minimum_maxdist=minimum_maxdist,
                          bad_bins=bad_bins)
  filters <- list()
  filters$border_size <- border_size
  filters$overhang <- overhang
  filters$AMP_THR <- AMP_THR
  filters$LOW_THR <- LOW_THR
  filters$MIN_WIDTH <- MIN_WIDTH
  filters$MIN_SEGMEAN_DIFF <- MIN_SEGMEAN_DIFF
  filters$min.gapwidth <- min.gapwidth
  filters$maxgap <- maxgap
  filters$minimum_maxdist=minimum_maxdist
  filters$minimum_count=minimum_count
  filters$edge <- edge
  filters
}

overlapsGermline <- function(object, germline, overhang=5e3){
  filter_levels <- c("no_germline", "part_germline", "fully_germline")
  start(germline) <- pmax(start(germline)-overhang, 1)
  sl <- seqlengths(germline)[as.character(seqnames(germline))]
  end(germline) <- pmin(end(germline)+overhang, sl)
  germline <- reduce(germline)
  any_overlap <- overlapsAny(object, germline)
  filtered <- ifelse(any_overlap, "part_germline", "no_germline")
  filtered[any_overlap] <- "part_germline"
  complete_overlap <- findOverlaps(object, germline, type="within")
  filtered[unique(queryHits(complete_overlap))] <- "fully_germline"
  filtered <- factor(filtered, levels=filter_levels)
  filtered
}

standardizeGRangesMetadata <- function(granges){
  if(length(granges) == 0) return(granges)
  mns <- granges$seg.mean
  is_amplicon <- granges$is_amplicon
  granges <- reduce(granges)
  granges$seg.mean <- mns
  granges$is_amplicon <- is_amplicon
  names(granges) <- ampliconNames(granges)

  granges$hgnc <- as.character(NA)
  granges$driver <- as.character(NA)
  granges$biol_sign <- as.character(NA)
  granges$groups <- as.factor(NA)
  granges
}

isNotAmplicon <- function(object) !isAmplicon(object)
isGap <- function(object) width(object)==999 & is.na(object$seg.mean)

gapBetweenLowCopyRanges <- function(g){
  index_lowcopy <- which(isNotAmplicon(g))
  index_gap1kb <- which(isGap(g))
  flanked_by_lowcopy <- (index_gap1kb - 1) %in% index_lowcopy &
    (index_gap1kb + 1) %in% index_lowcopy
  gap1kb <- g[index_gap1kb]
  gap1kb[flanked_by_lowcopy]
}

startAmpliconBoundary <- function(ranges, is_gap=TRUE, extend=2e3){
  if(!is_gap){
    candidate_amplicon <- ranges$seg.mean > 0.8 & !is.na(ranges$seg.mean)
    if(any(!candidate_amplicon)){
      start(ranges)[!candidate_amplicon] <- end(ranges)[!candidate_amplicon]
    }
  }
  is_small <- width(ranges) < extend
  if(any(is_small)){
    small_amplicons <- ranges[is_small]
    if(all(is_small)) return(small_amplicons)
  }
  amplicons <- ranges[!is_small]
  cnv_begin_edge <- amplicons
  start(cnv_begin_edge) <- start(amplicons)
  end(cnv_begin_edge) <- start(amplicons)+extend
  if(any(is_small)){
    cnv_begin_edge <- sort(c(cnv_begin_edge, small_amplicons))
  }
  cnv_begin_edge
}

endAmpliconBoundary <- function(ranges, is_gap=FALSE, extend=2e3){
  if(!is_gap){
    candidate_amplicon <- ranges$seg.mean > 0.8 & !is.na(ranges$seg.mean)
    if(any(!candidate_amplicon)){
      start(ranges)[!candidate_amplicon]<- end(ranges)[!candidate_amplicon]
    }
  }
  is_small <- width(ranges) < extend
  if(any(is_small)){
    small_amplicons <- ranges[is_small]
    if(all(is_small)) return(small_amplicons)
  }
  amplicons <- ranges[!is_small]
  if(length(amplicons) == 0) return(amplicons)
  cnv_end_edge <- amplicons
  start(cnv_end_edge) <- end(amplicons)-extend
  end(cnv_end_edge) <- end(amplicons)
  if(any(is_small)){
    cnv_end_edge <- sort(c(cnv_end_edge, small_amplicons))
  }
  cnv_end_edge
}

borderSize <- function(object) object@border_size

germlineCNV <- function(object) object@germline_cnv

outliers <- function(object) object@outliers

ampliconQueryRanges <- function(object, min.fc=1){
  g <- ranges(object)
  gap1kb <- gapBetweenLowCopyRanges(g)
  starts <- startAmpliconBoundary(g, extend=borderSize(object))
  ends <- endAmpliconBoundary(g, extend=borderSize(object))
  g <- c(starts, ends)
  names(g) <- ampliconNames(g)
  g <- g[!duplicated(names(g))]
  g <- g[width(g) > 1]
  g <- g[!is.na(g$seg.mean)]
  g <- g[g$seg.mean > min.fc | isGap(g)]
  ## drop 999bp intervals between low copy ranges
  g <- g[!names(g) %in% names(gap1kb)]
  assembly_gaps <- assemblyGaps(object)
  if(length(assembly_gaps)==0) {
    queryRanges(object) <- g
    return(object)
  }
  ## drop ranges contained within assembly gaps
  assembly_gaps <- reduce(assembly_gaps)
  o <- findOverlaps(g, assembly_gaps, type="within")
  if(length(o) > 0){
    g <- g[-unique(queryHits(o))]
    if(length(g)==0) warning("All queryRanges are within assembly gaps")
  }
  ##  ## drop ranges contained within germline CNV / outliers gaps
  germline <- sort(reduce(c(germlineCNV(object), outliers(object))))
  o <- findOverlaps(g, germline, type="within")
  if(length(o) > 0){
    g <- g[-unique(queryHits(o))]
    if(length(g)==0) warning("All queryRanges are within germline")
  }
  queryRanges(object) <- g
  object
}

setMethod("combine", signature(x="GRanges", y="GRanges"),
          function(x, y,...){
            col.names <- intersect(colnames(mcols(x)), colnames(mcols(y)))
            ##col.names <- union(colnames(mcols(x)), colnames(mcols(y)))
            mcols(x) <- mcols(x)[, col.names]
            mcols(y) <- mcols(y)[, col.names]
            colnames(mcols(x)) <- colnames(mcols(y)) <- col.names
            xy <- c(x, y)
            xy
          })

.add_mcols <- function(g){
  g$seg.mean <- NA
  g$is_amplicon <- FALSE
  g$hgnc <- as.character(NA)
  g$driver <- as.character(NA)
  g$biol_sign <- as.character(NA)
  g$groups <- as.factor(NA)
  g
}

.empty_filters <- function(assembly_gaps=GRanges(),
                           centromeres=GRanges(),
                           germline_cnv=GRanges(),
                           outliers=GRanges()){
  list(assembly_gaps=assembly_gaps,
       centromeres=centromeres,
       germline_cnv=germline_cnv,
       outliers=outliers)
}

#' Constructor for AmpliconGraph
#'
#' Constructor for \code{AmpliconGraph}.
#'
#' @param ranges a \code{GRanges} of the amplicons
#' 
#' @param border_size used to construct a query for additional
#'   amplicons neighboring a focal amplicon. TODO: more detail needed.
#' 
#' @param assembly_gaps a \code{GRanges} object of the assembly gaps
#' 
#' @param centromeres a \code{GRanges} object of the centromeres
#' 
#' @param germline_cnv a \code{GRanges} object of germline CNVs
#' 
#' @param outliers a \code{GRanges} object of germline outliers
#' 
#' @param overhang a length-one numeric vector
#' @rdname AmpliconGraph-constructor
#' @export
AmpliconGraph <- function(ranges=GRanges(),
                          filters,
                          params=ampliconParams(border_size=150e3, overhang=5e3)){
  ##                          border_size=150e3,
  ##                        assembly_gaps=GRanges(),
  ##                        centromeres=GRanges(),
  ##                        germline_cnv=GRanges(),
  ##                        outliers=GRanges(),
  ##                        overhang=5e3){
  if(missing(filters)){
    filters <- .empty_filters()
  }
  assembly_gaps <- filters[["assembly_gaps"]]
  centromeres <- filters[["centromeres"]]
  germline_cnv <- filters[["germline_cnv"]]
  outliers <- filters[["outliers"]]
  border_size <- params[["border_size"]]
  overhang <- params[["overhang"]]
  if(missing(ranges)) {
    ag <- new("AmpliconGraph")
    names(ranges(ag)) <- character()
    return(ag)
  }
  ranges <- standardizeGRangesMetadata(ranges)
  assembly_gaps <- reduce(sort(combine(assembly_gaps, centromeres)))
  ##
  ## Annotate GRanges
  ##
  ##  gaps
  ##
  g <- gaps0(ranges)
  if(length(g) > 0){
    g <- .add_mcols(g)
    ranges <- sort(c(ranges, g))
  }
  ##
  germ <- reduce(c(germline_cnv, outliers))
  ranges$overlaps_germline <- overlapsGermline(ranges, germ,
                                               overhang=overhang)
  names(ranges) <- ampliconNames(ranges)
  amps <- ranges[isAmplicon(ranges) & isSomatic(ranges)]
  tmp <- new("AmpliconGraph",
             graph=graphNEL(nodes=names(amps)),
             ranges=ranges,
             border_size=border_size,
             assembly_gaps=assembly_gaps,
             germline_cnv=germline_cnv,
             outliers=outliers)
  ## initialize query ranges
  ampliconQueryRanges(tmp)
}

trimRangesOverlappingCentromere <- function(object, centromeres){
  if(numEdges(graph(object)) > 0) {
    stop("currently this method does not preserve edges")
  }
  rg <- ranges(object)
  sl <- seqlevels(rg)
  si <- seqinfo(rg)
  spanned_by_gap <- unique(queryHits(findOverlaps(rg, centromeres, type="within")))
  if(length(spanned_by_gap) > 0){
    rg <- rg[-spanned_by_gap]
  }
  overlaps_centromere <- overlapsAny(rg, centromeres)
  if(any(overlaps_centromere)){
    rgs <- rg[overlaps_centromere]
    trimmed <- setdiff(rgs, centromeres)
    j <- subjectHits(findOverlaps(trimmed, rgs))
    ## only keep one side
    tmp <- sapply(split(trimmed, j), function(g) g[which.max(width(g))])
    trimmed <- unlist(GRangesList(tmp))
    j <- subjectHits(findOverlaps(trimmed, rgs))
    trimmed$seg.mean <- rgs$seg.mean[j]
    trimmed$is_amplicon <- rgs$is_amplicon[j]
    trimmed$hgnc <- rgs$hgnc
    trimmed$driver <- rgs$driver
    trimmed$biol_sign <- rgs$biol_sign
    trimmed$groups <- rgs$groups
    trimmed$overlaps_germline <- rgs$overlaps_germline
    ## drop regions that were gaps to begin with
    trimmed <- trimmed[!is.na(trimmed$seg.mean)]
    rg <- filterBy(rg, rgs)
    rg <- sort(c(rg, trimmed))
    names(rg) <- ampliconNames(rg)
    rg <- keepSeqlevels(rg, sl, pruning.mode="coarse")
    seqinfo(rg) <- si
    ranges(object) <- rg
  }
  ag <- graphNEL(nodes=names(ampliconRanges(object)))
  graph(object) <- ag
  object
}

#' Merge adjacent amplicons that have a similar segment mean
#'
#' Merges adjacent amplicons (<= 3kb separation) that have similar
#' segment means.
#'
#' @keywords internal
#' @export
#' @return a \code{GRanges} object of the merged segments
#'
#' @param object a \code{GRanges} object
#'
#' @param thr a length-one numeric vector. If the difference of the average log
#'   ratio for two segments is less than this number, combine the segments. The
#'   segment mean of the merged GRanges will be the weighted average, where the
#'   weights are the segment widths.
#' 
#'
#' @param MIN_WIDTH minimum size of amplicon (default 2000)
#'
#' @details Merging is performed by the \code{reduce} operation.
#' 
#'
joinNearGRanges <- function(object, params){
  ##thr=0.05, MIN_WIDTH=2e3, min.gapwidth=3000){
  thr <- params[["MIN_SEGMEAN_DIFF"]]
  MIN_WIDTH <- params[["MIN_WIDTH"]]
  min.gapwidth <- params[["min.gapwidth"]]
  k <- which(!is.na(object$seg.mean) & width(object) > MIN_WIDTH)
  if(length(k) == 0) return(object)
  g <- object[k]
  means <- g$seg.mean
  d <- diff(means)
  if(all(d > thr)) return(object)
  regions <- paste0("region", cumsum(c(0, abs(d) > thr)))
  regions <- factor(regions, levels=unique(regions))
  indexlist <- split(seq_along(regions), regions)
  ##
  ## list of regions.  Each element is a vector of regions that can be
  ## joined
  ##
  glist <- vector("list", length(indexlist))
  for(i in seq_along(indexlist)){
    j <- indexlist[[i]]
    if(length(j)==1) next()
    ng <- reduce(g[j], min.gapwidth=min.gapwidth)
    ##
    ## TODO: combining the metadata.  Here, we're just taking the first row
    ##
    mc <- mcols(g)[j[1], , drop=FALSE]
    mc$seg.mean <- sum(width(g[j])*g$seg.mean[j])/sum(width(g[j]))
    mcols(ng) <- mc
    glist[[i]] <- ng
  }
  gnew <- unlist(GRangesList(unlist(glist)))
  ##  Remove regions in the original object that overlap with the
  ## joined regions (g contains all the original intervals that were
  ## subsequently joined)
  object2 <- filterBy(object, gnew)
  object3 <- c(object2, gnew)
  object3$seg.mean <- round(object3$seg.mean, 2)
  object3 <- sort(object3)
  object3
}

flankingRanges <- function(object, shift=2){
  rgs <- ranges(object)
  amplicons <- ampliconRanges(object)
  index <- match(names(amplicons), names(rgs))
  index2 <- index - shift  
  index2 <- index2[index2 > 0]
  index3 <- index + shift
  is_right_flank <- index3 < length(rgs)
  if(any(!is_right_flank)){
    ## REFACTOR.  No right flank exists.  If we exclude, linkedDuplicatedRanges will fail
    ## assign the amplicon index
    index3[!is_right_flank] <- index[!is_right_flank]
  }  
  left_flank <- rgs[index2]
  right_flank <- rgs[index3]
  left <- list(first=left_flank, last=amplicons)
  right <- list(first=amplicons, last=right_flank)
  list(left=left, right=right)
}

flankingDuplications <- function(object, minimum_foldchange=1){
  flanks <- flankingRanges(object,  shift=2)
  is_dupL <- isDuplication(flanks$left[[1]], minimum_foldchange)
  is_dupR <- isDuplication(flanks$right[[2]], minimum_foldchange)
  flanks$left <- lapply(flanks$left, "[", is_dupL)
  flanks$right <- lapply(flanks$right, "[", is_dupR)
  flanks
}

minimumBasepairCoverage <- function(rp, regions){
  dj <- disjoin(rp)
  dj <- dj[overlapsAny(dj, regions)]
  cnt <- countOverlaps(dj, rp)
  h <- findOverlaps(regions, dj)
  cntlist <- split(cnt[subjectHits(h)], queryHits(h))
  min_cnt <- sapply(cntlist, min)
  min_cnt
}

linkedDuplicatedRanges <- function(object, rpsegs,
                                   flanking_duplications,
                                   params){
  minimum_count <- params[["min_count"]]
  flank <- flanking_duplications
  lengths <- unlist(lapply(flank, elementNROWS))
  ##if(any(lengths != 1)) stop("Flanking regions must be length-one GRanges")
  gapsLeft <- GRanges(seqnames(flank[["left"]]$first),
                      IRanges(end(flank[["left"]]$first)+1,
                              start(flank[["left"]]$last)))
  gapsRight <- GRanges(seqnames(flank[["right"]]$first),
                       IRanges(end(flank[["right"]]$first)+1,
                               start(flank[["right"]]$last)))
  cntsLeft <- minimumBasepairCoverage(rpsegs, gapsLeft)
  cntsRight <- minimumBasepairCoverage(rpsegs, gapsRight)
  left <- lapply(flank$left, "[", which(cntsLeft >= minimum_count))
  right <- lapply(flank$right, "[", which(cntsRight >= minimum_count))
  list(left=left, right=right)##flanking_duplications[cnts >= minimum_count]
}

addFlanks <- function(object, dup_granges){
  ## left-of-amplicon flank
  flank <- dup_granges$left[[1]]
  if(length(flank) > 0){
    amp <- dup_granges$left[[2]]
    edges <- paste(names(flank), names(amp), sep="-")
    index <- match(names(flank), names(ranges(object)))
    new_nodes <- names(ranges(object))[index]
    object <- addNode(new_nodes, object, edges)
  }
  ## right-of-amplicon flank
  flank <- dup_granges$right[[2]]
  if(length(flank) > 0){
    amp <- dup_granges$right[[1]]
    edgR <- paste(names(amp), names(flank), sep="-")
    index <- match(names(flank), names(ranges(object)))
    new_nodes <- names(ranges(object))[index]
    object <- addNode(new_nodes, object, edgR)
  }
  object
}

#' Add focal amplicons that flank an amplicon seed 
#'
#' Amplicons selected to seed a graph have a fold change of at least log2(2.75),
#' or nearly 3-fold the diploid genome. Lower-copy amplicons flanking the seed
#' amplicons are added to the graph if they have a fold change of at least
#' log2(1.75), or nearly 2-fold the diploid genome. To assess whether a low-copy
#' amplicon is 'flanking' an amplicon seed, all read pairs that aligned to the
#' amplicon seeds (see \code{get_readpairs}) are passed to this function in
#' object \code{rp}. The read pairs are converted to segments by
#' \code{readPairsAsSegments} and an assessment of whether flanking duplications
#' are linked is given by \code{linkedDuplicatedRanges}. In particular, this
#' function checks whether there are at least \code{minimum_count} fragments
#' that connect the amplicon seeds to the flanking regions. Finally, the flanks
#' are added to the graph by the \code{addFlanks} function.
#'
#' REFACTORING: overly complicated with many abstractions.  Needs an
#' example to step through
#'
#'
#' @param object a \code{AmpliconGraph}
#'
#' @param rp a \code{GAlignmentPairs} object
#'
#' @param params a list of parameters for amplicon analyses
#' @seealso \code{\link{ampliconParams}}}
#' @examples
#'   path <- file.path(system.file(package="svcnvs"), "tests/testthat")
#'   ag <- readRDS(file.path(path, "addFocalDups.ffab104.rds"))
#'   params <- ampliconParams()
#'   ag <- addFocalDupsFlankingAmplicon(ag, rp, params)
#' @export
addFocalDupsFlankingAmplicon <- function(object, rp, params){
  if(totalWidth(queryRanges(object))==0) return(object)
  flanks <- flankingDuplications(object,
                                 minimum_foldchange=params[["LOW_THR"]])
  rpsegs <- readPairsAsSegments(rp)
  dup_gr <- linkedDuplicatedRanges(object, rpsegs, flanks, params)
  object <- addFlanks(object, dup_gr)
  object
}

numberAnchored <- function(object, read_pairs){
  seed <- ampliconRanges(object)
  R1_anchored <- overlapsAny(first(read_pairs), seed)
  R2_anchored <- overlapsAny(last(read_pairs), seed)
  number_anchored <- R1_anchored+R2_anchored
  zero_anchor <- read_pairs[number_anchored == 0]
  single_anchor <- read_pairs[number_anchored == 1]
  both_anchored <- read_pairs[number_anchored == 2]
  AnchoredReadPairs(zero_end=zero_anchor,
                    single_end=single_anchor,
                    both_ends=both_anchored)
}

edgeStats <- function(edges, param){
  freq <- table(names(edges))
  freq <- freq[freq >= param$freq]
  if(length(freq)==0) return(data.frame(freq=integer(), maxdist=integer(), bad_bin=logical()))
  stats <- data.frame(freq=as.integer(freq))
  rownames(stats) <- names(freq)
  e2 <- edges[names(edges) %in% rownames(stats)]
  index <- split(seq_len(length(e2)), names(e2))
  ##i <- NULL
  results <- vector("list", length(index))
  for(j in seq_along(index)){
    i <- index[[j]]
    results[[j]] <- start(first(e2))[i]
  }
  maxdist <- sapply(lapply(results, range), diff)
  ##maxdist <- sapply(lapply(foreach(i = index) %do% start(first(e2))[i], range), diff)
  names(maxdist) <- names(index)
  stats$maxdist <- maxdist[rownames(stats)]
  stats$bad_bin <- rep(FALSE, nrow(stats))
  bad_bins <- param$bad_bins
  if(length(bad_bins) == 0 ) return(stats)
  hits1 <- findOverlaps(first(e2), bad_bins[[1]], select="first")
  hits2 <- findOverlaps(last(e2), bad_bins[[2]], select="first")
  if(any(hits1 == hits2, na.rm=TRUE)){
    ## matches a flagged bin pair
    stop("variable flag not defined")
    ##id <- names(e2)[flag]
    ##stats$bad_bin[id] <- TRUE
  }
  ## check the reverse
  hits1 <- findOverlaps(first(e2), bad_bins[[2]], select="first")
  hits2 <- findOverlaps(last(e2), bad_bins[[1]], select="first")
  if(any(hits1 == hits2, na.rm=TRUE)){
    stop("variable flag not defined")
##    id <- names(e2)[flag]
##    stats$bad_bin[id] <- TRUE
  }
  stats
}

edgeFilters <- function(edges, param=FilterEdgeParam()){
  stats <- edgeStats(edges, param)
  valid_edges <- rownames(stats)[stats$maxdist >= param$minimum_maxdist & !stats$bad_bin]
  valid_edges
}

singleAnchorHits <- function(object, read_pair){
  if(length(read_pair) < 1) return(NULL)
  a <- ampliconRanges(object)
  hitfirst <- findOverlaps(first(read_pair), a, select="all", maxgap=1e3)
  g <- nonampliconRanges(object)
  hitlast <- findOverlaps(last(read_pair), g, select="all", maxgap=1e3)
  hitlist <- hitsWithSameQuery(list(hitfirst, hitlast))
  if(!identicalReadQueries(hitlist)){
    hitlist <- dropDuplicatedQueryHits(hitlist)
  }
  hitlist <- hitsWithDifferentTargets(hitlist)
  edges <- orderEdgesByIndexOfOverlap(hitlist, ranges1=a, ranges2=g)
  supporting_read_pair <- read_pair[queryHits(hitlist[[1]])]
  names(supporting_read_pair) <- edges
  supporting_read_pair
}

addDuplications <- function(object, granges, edges){
  if(length(granges) > 0){
    index <- match(unique(names(granges)), names(ranges(object)))
    new_nodes <- names(ranges(object))[index]
    object <- addNode(new_nodes, object, edges)
  }
  object
}

linkFocalDups <- function(object, rp, params){
  LOW_THR <- params[["LOW_THR"]]
  edgeParam <- params[["edge"]]
  if(totalWidth(queryRanges(object))==0) return(object)
  anchors <- numberAnchored(object, rp)
  se_anchors <- singleEndAnchors(anchors)
  e <- singleAnchorHits(object, se_anchors)
  if(is.null(e)) return(object)
  keep <- edgeFilters(e, param=edgeParam)
  e <- e[names(e) %in% keep]
  if(length(e)==0) return(object)
  candidates <- ranges(object)[node2(names(e))]
  is_dup <- isDuplication(candidates, minimum_foldchange=LOW_THR)
  e <- e[is_dup]
  candidates <- candidates[is_dup]
  object <- addDuplications(object, candidates, unique(names(e)))
  object
}

orderEdgesByIndexOfOverlap <- function(hitlist, ranges1, ranges2){
  if(missing(ranges2)) ranges2 <- ranges1
  j <- subjectHits(hitlist[[1]])
  k <- subjectHits(hitlist[[2]])
  jk <- cbind(pmin(j, k), pmax(j,k))
  edges <- paste(names(ranges1)[jk[,1]],
                 names(ranges2)[jk[,2]],
                 sep="-")
  edges
}

hitsWithSameQuery <- function(hitlist){
  qhits <- intersect(queryHits(hitlist[[1]]), queryHits(hitlist[[2]]))
  hitlist <- lapply(hitlist, function(x, index) x[queryHits(x) %in% index], index=qhits)
  hitlist
}

hitsWithDifferentTargets <- function(hitlist){
  j <- subjectHits(hitlist[[1]])
  k <- subjectHits(hitlist[[2]])
  notequal <- j != k
  hitlist <- lapply(hitlist, "[", i=notequal)
  hitlist
}

identicalReadQueries <- function(hitlist)
  identical(queryHits(hitlist[[1]]), queryHits(hitlist[[2]]))

dropDuplicatedQueryHits <- function(hitlist){
  qh1 <- queryHits(hitlist[[1]])
  qh2 <- queryHits(hitlist[[2]])
  d1 <- duplicated(qh1)
  d2 <- duplicated(qh2)
  hitlist[[1]] <- hitlist[[1]][!d1, ]
  hitlist[[2]] <- hitlist[[2]][!d2, ]
  hitlist
}

twoAnchorHits <- function(object, read_pair){
  if(length(read_pair) < 1) return(NULL)
  a <- ampliconRanges(object)
  hitfirst <- findOverlaps(first(read_pair), a, select="all", maxgap=1e3)
  hitlast <- findOverlaps(last(read_pair), a, select="all", maxgap=1e3)
  hitlist <- hitsWithSameQuery(list(hitfirst, hitlast))
  if(!identicalReadQueries(hitlist)){
    hitlist <- dropDuplicatedQueryHits(hitlist)
  }
  hitlist <- hitsWithDifferentTargets(hitlist)
  edges <- orderEdgesByIndexOfOverlap(hitlist, ranges1=a)
  supporting_read_pair <- read_pair[queryHits(hitlist[[1]])]
  names(supporting_read_pair) <- edges
  supporting_read_pair
}

linkAmplicons <- function(object, rp, edgeParam=FilterEdgeParam()){
  rp <- rp[overlapsAny(rp, ampliconRanges(object), maxgap=1e3)]
  if(numNodes(object)<2) return(object)
  anchors_rp <- numberAnchored(object, rp)
  e <- twoAnchorHits(object, bothEndAnchors(anchors_rp))
  if(length(e)==0) return(object)  ##there are no improper pairs linking amplicons
  e <- e[names(e) %in% edgeFilters(e, param=edgeParam)]
  if(length(e)==0) return(object)
  tabe <- table(names(e))
  keep <- names(tabe)[tabe >= 5]
  from <-  node1(keep)
  to <- node2(keep)
  existing <- edges(object)[from] ## is 'to' in any of the existing edges
  edge_exists <- mapply(function(to, existing) to %in% existing,
                        to=to, existing=existing)
  if(!all(edge_exists)){
    graph(object) <- addEdge(node1(keep[!edge_exists]),
                         node2(keep[!edge_exists]), graph(object))
  }
  object
}

linkNearAmplicons <- function(object, maxgap=500e3){
  hits <- findOverlaps(ampliconRanges(object), maxgap=maxgap)
  hits <- hits[!isSelfHit(hits)]
  hits <- hits[!isRedundantHit(hits)]
  if(length(hits)==0) return(object)
  new_edges <- paste(names(ampliconRanges(object))[queryHits(hits)],
                     names(ampliconRanges(object))[subjectHits(hits)],
                     sep="-")
  from <-  node1(new_edges)
  to <- node2(new_edges)
  existing <- edges(object)[from] ## is 'to' in any of the existing edges
  edge_exists <- mapply(function(to, existing) to %in% existing,
                        to=to, existing=existing)
  new_edges <- new_edges[!edge_exists]
  if(length(new_edges)==0) return(object)
  graph(object) <- addEdge(node1(new_edges),
                           node2(new_edges), graph(object))
  object
}

germline <- function(object){
  reduce(c(germlineCNV(object), outliers(object)))
}

#' Identify lower copy focal amplicons (possibly duplicatons) not already
#' annotated as amplicons and not contained within the germline filters
#'
#' Select all duplicated segments (see \code{isDuplication}) less than `maxgap`
#' kb from the seed amplicons of the graph that are not within a germline CNV or
#' outlier. Germline CNVs and outliers can be conveniently extracted from the
#' graph object as given in the examples.
#'
#' @export
#'
#' @param object A \code{AmpliconGraph}
#' @param params a list of parameters such as that generated by \code{ampliconParams}.
#' @return a reduced \code{GRanges} of lower-copy amplicons
#' @seealso \code{\link{ampliconParams}}}
#' @examples
#'   ## load a previously saved AmpliconGraph from a regression unit test
#'   path <- system.file("testthat", package="svcnvs")
#'   ag <- readRDS(file.path(path, "addFocalDups.ffab104.rds"))
#'   params <- ampliconParams()
#'   focalAmpliconDupRanges(ag, params)
#' @export
focalAmpliconDupRanges <- function(object, params){
  LOW_THR <- params[["LOW_THR"]]
  MAX_SIZE <- params[["maxgap"]]
  g <- ampliconRanges(object)
  is_dup <- isDuplication(ranges(object), minimum_foldchange=LOW_THR)
  ## remove any lower copy ranges that are very large
  dup_g <- ranges(object)[is_dup]
  dup_g <- dup_g[width(dup_g) < MAX_SIZE]
  g <- sort(c(dup_g, g))
  ##
  ## TODO: add expansion to params
  ##
  if(length(g) > 0){
    g <- expandGRanges(g, 1e3)
  }
  ## we do not want to link germline events
  germ <- germline(object)
  ##
  ## TODO: add expansion to params
  ##
  germ <- expandGRanges(germ, 10e3)
  germ <- reduce(germ)
  g <- filterBy(g, germ, type="within")
  reduce(g)
}

removeNode2 <- function (node, object) {
  ##object <- clearNode(node, object)
  nN <- nodes(object)
  node <- node[node %in% nN]
  if(length(node) == 0) return(object)
  wh <- match(node, nN)
  gN <- nN[-wh]
  nE <- object@edgeL
  nE <- nE[gN]
  object@edgeL <- nE
  identical(names(nE), gN)
  ##nE <- object@edgeL[-wh]
  nE2 <- lapply(nE, function(el) {
    ## Nodes that never linked to anything to begin with
    if(length(el$edges) == 0) return(el)
    oldN <- nN[el$edges]
    result <- match(oldN, gN)
    result <- as.integer(result[!is.na(result)])
    el$edges <- result
    el
  })
  object@nodes <- gN
  object@edgeL <- nE2
  object
}

filterSmallAmplicons <- function(object, MIN_SIZE=5e3){
  ar <- reduce(ampliconRanges(object), min.gapwidth=5e3)
  ar <- ar[width(ar) < 5e3]
  if(length(ar) > 0){
    j <- subjectHits(findOverlaps(ar, ranges(object)))
    gN <- removeNode2(names(ranges(object))[j], graph(object))
    graph(object) <- gN
  }
  object
}

updateRangesMetadata <- function(object, granges){
  index <- match(names(granges), names(ranges(object)))
  rg <- ranges(object)[-index]
  rg2 <- sort(c(rg, granges))
  ranges(object) <- rg2
  object
}

#' Assign a group id for nodes that are linked
#'
#' Takes a list of amplicons (GRanges object) and an undirected graph
#' representation of the amplicons
#'
#' @export
#' @param gN a \code{AmpliconGraph} object
groupAmplicons <- function(gN){
  V <- nodes(gN)
  if(numEdges(gN) > 0){
    kgroups <- kCliques(gN)
    groups <- kgroups[[length(kgroups)]]
    grouped.amplicon <- rep(NA, length(V))
    for(k in seq_along(groups)){
      grouped.amplicon[match(groups[[k]], V)] <- k
    }
    if(any(is.na(grouped.amplicon))){
      grouped.amplicon[is.na(grouped.amplicon)] <- max(grouped.amplicon, na.rm=TRUE) +
        seq_len(sum(is.na(grouped.amplicon)))
    }
    grouped.amplicons <- as.integer(factor(grouped.amplicon))
  } else grouped.amplicon <- seq_len(length(V))
  names(grouped.amplicon) <- V
  factor(grouped.amplicon)
}

setAmpliconGroups <- function(object){
  groups <- groupAmplicons(graph(object))
  groups <- setNames(paste0("gr", groups), names(groups))
  ar <- ampliconRanges(object)
  groups <- groups[names(ar)]
  ar$groups <- groups
  object <- updateRangesMetadata(object, ar)
  object
}

getGenes <- function(object, transcripts){
  .hgnc <- setNames(rep(NA, length(object)), names(object))
  hits <- findOverlaps(transcripts, object, type="within")
  ##genes <- hgnc(transcripts)[queryHits(hits)]
  genes <- transcripts$gene_name[queryHits(hits)]
  amp <- names(object)[subjectHits(hits)]
  geneL <- lapply(split(genes, amp), function(x) paste(unique(x), collapse=", "))
  genes <- unlist(geneL)
  .hgnc[names(genes)] <- genes
  object$hgnc <- as.character(.hgnc)
  object
}

setGenes <- function(object, transcripts){
  ar <- getGenes(ampliconRanges(object), transcripts)
  object <- updateRangesMetadata(object, ar)
  object
}

## Two levels of significance for genes:
##    - clinical significance (synonymous with driver)
##
##    - biologically significant: perhaps biologically significant but unkown
##      clinical significance. This set of all clinically significant
##      genes is a subset
driver_genes <- function(tx, clin_sign=FALSE){
  if(clin_sign){
    return(tx$gene_name[tx$cancer_connection])
  }
  tx$gene_name [ tx$biol_sign ]
}

getDrivers <- function(object, transcripts, clin_sign=TRUE){
  known_drivers <- unique(driver_genes(transcripts, clin_sign=clin_sign))
  genes <- object$hgnc
  gene_list <- split(genes, object$groups)
  gene_list <- lapply(gene_list, function(x){
    x <- x[!is.na(x)]
    if(length(x)==0) return(NA)
    xx <- unlist(strsplit(x, ", "))
  })
  driver_in_grouped_amplicon <- lapply(gene_list, function(x){
    x <- paste(x[x %in% known_drivers], collapse=", ")
    if(x == "") x <- NA
    x
  })
  driver_in_grouped_amplicon <- unlist(driver_in_grouped_amplicon)
  if(all(is.na(driver_in_grouped_amplicon))) return(object)
  drv <- driver_in_grouped_amplicon[!is.na(driver_in_grouped_amplicon)]
  for(i in seq_along(drv)){
    driver_group <- names(drv)[i]
    if(clin_sign){
      object$driver[object$groups == driver_group] <- drv[i]
    } else {
      object$biol_sign[ object$groups == driver_group ] <- drv[i]
    }
  }
  if(clin_sign){
    object$driver <- as.character(object$driver)
  } else {
    object$biol_sign <- as.character(object$biol_sign)
  }
  object
}

setDrivers <- function(object, transcripts, clin_sign=TRUE){
  ar <- ampliconRanges(object)
  ar <- getDrivers(ar, transcripts, clin_sign=clin_sign)
  object <- updateRangesMetadata(object, ar)
  object
}

makeAGraph <- function(segs, af, params){
  segs$is_amplicon <- segs$seg.mean > params$AMP_THR
  ag <- AmpliconGraph(ranges=segs,
                      filters=af,
                      params=params)
  ##border_size=af[["border_size"]],
  ##                    assembly_gaps=af[["assembly_gaps"]],
  ##                    centromeres=af[["centromeres"]],
  ##                    germline_cnv=af[["germline_cnv"]],
  ##                    outliers=af[["outliers"]],
  ##                    overhang=af[["overhang"]])
  if (numNodes (ag) == 0) return (ag)
  ag <- trimRangesOverlappingCentromere (ag, af[["centromeres"]])
  ag
}

#' Construct an AmpliconGraph from a BamViews object
#'
#' This function constructs an \code{AmpliconGraph} from a
#' \code{BamViews} object for a single sample.  By default, the seeds
#' of the graph are focal amplicons with fold-change of nearly 3
#' relative to the diploid genome (log2(2.75)).  The threshold of
#' seeding amplicons can be adjusted by the \code{ampliconParams}
#' function. After seeding the graph with high-copy focal amplicons,
#' both neighboring (flanking) and distant low-copy focal amplicons
#' are added to the graph object.  Next, improperly paired reads in
#' which both the first and last read align to any queryRange of the
#' graph object are parsed from the bam file ( this function will
#' throw an error if not all files in \code{bamPaths} exist).  If 5 or
#' more improperly paired reads bridge a node to another node, these
#' amplicons are grouped.  Further, if a low-copy amplicon is bridged
#' to an existing node, the low-copy amplicon will become a node in
#' the graph. Amplicon groups are defined by the edges between nodes,
#' where the edges represent improperly paired reads that support a
#' junction between two amplicons.
#'
#' REFACTORING: needs an example to step through.  Perhaps an initial
#' graph object and then keep updating the graph object and associated
#' visualization with each step.  Each of the low level functions in
#' sv_deletions should be exported to more fully document this
#' procedure.
#'
#' @seealso See \code{ampliconParams} for default parameters. The
#'   wrapper \code{\link{sv_amplicon_exp}} constructs and saves an
#'   \code{\linkS4class{AmpliconGraph}} for each sample in an
#'   experiment.
#'
#' @examples
#'   library(svovarian)
#'   library(svfilters.hg19)
#'   id <- "CGOV2T"
#'   data(lymph_ids)
#'   data(germline_filters)
#'   data(transcripts)
#'   params <- ampliconParams()
#'   bviews <- readRDS(file.path(dp[1], "bviews_hg19.rds"))
#'   bviews <- bviews[, id]
#'   grl <- readRDS(file.path(dp["segment"], "grl_hg19.rds"))
#'   ## Requires bam file
#'   ag <- sv_amplicons(bviews[, id], grl[[id]], germline_filters, params)
#'
#' @export
#' @param bview a \code{BamViews} object
#'
#' @param segs a \code{GRanges} object of segments with log2
#'   fold-changes consistent with amplification
#'
#' @param amplicon_filters a \code{list} of filters
#' @param params a list of parameters for the amplicon analysis
#' @param transcripts a \code{GRanges} object of transcripts
sv_amplicons <- function(bview, segs, amplicon_filters, params, transcripts){
  ag <- makeAGraph(segs, amplicon_filters, params)
  tmp <- joinNearGRanges(ranges(ag), params)
  names(tmp) <- ampliconNames(tmp)
  ## the names of the nodes no longer correspond to the range names
  ## stopifnot(nodes(ag) %in% names(tmp)), and so
  ## setAmpliconGroups fails
  ranges(ag) <- tmp
  REMOTE <- file.exists(bamPaths(bview))
  if(!REMOTE) stop ("Path to BAM files is invalid")
  LOW_THR <- params[["LOW_THR"]]
  ##
  ## REFACTOR: could this step use the saved improper read pairs
  ##
  rp <- get_readpairs(ag, bamPaths(bview))
  ag <- addFocalDupsFlankingAmplicon(ag, rp, params)
  queryRanges(ag) <- focalAmpliconDupRanges(ag, params)
  irp <- get_improper_readpairs(ag, bamPaths(bview))
  ##
  ## At this point, focal duplications added to the graph have not
  ## been linked to any of the seeds
  ##
  ##paired_bin_filter <- af[["paired_bin_filter"]]
  ##param <- FilterEdgeParam(minimum_maxdist=50, bad_bins=paired_bin_filter)
  ##param <- FilterEdgeParam(minimum_maxdist=50, bad_bins=GRanges())
  edge.p <- params[["edge"]]
  ag <- linkFocalDups(ag, irp, params)
  ag <- linkAmplicons(ag, irp, edgeParam=edge.p)
  ag <- linkNearAmplicons(ag, maxgap=params[["maxgap"]])
  ag <- filterSmallAmplicons (ag)
  ag <- setAmpliconGroups (ag)
  ag <- setGenes (ag, transcripts)
  ag <- setDrivers (ag, transcripts, clin_sign=TRUE)
  ag <- setDrivers (ag, transcripts, clin_sign=FALSE)
  ag
}



#' Identify somatic amplicons and grouped amplicons
#'
#' This function identifies focal, somatic amplifications.  As
#' amplicons are often grouped by the bridge-fusion-breakage cycles,
#' we represent the set of somatic amplicons identified for a sample
#' as a graph. The nodes of the graph are the amplicons and the edges
#' are the links between amplicons informed by paired reads. This
#' function is a wrapper for \code{sv_amplicons} that constructs an
#' \code{AmpliconGraph} for a given sample.  As with the other
#' *Experiment functions, \code{sv_amplicon_exp} saves the
#' \code{AmpliconGraph} computed for each sample as in intermediate
#' file for quick recall.  If the file exists, this function reads the
#' serialized R object from disk. To generate an \code{AmpliconGraph}
#' de novo, one must first delete the intermediate files (see
#' examples).
#'
#' @seealso See \code{\linkS4class{AmpliconGraph}} for methods
#'   associated with the class and \code{\link{sv_amplicons}} for
#'   construction of an \code{AmpliconGraph} for a single sample.
#' 
#'
#' @param dirs character-vector of file paths for storing intermediate
#'   files
#' @param bviews A \code{BamViews} object
#' @param grl A \code{GRangesList} of the segmented genomes (each
#'   element is the \code{GRanges} for a sample)
#' @param amplicon_filters A list of germline filters and parameters. See the \code{germline_filters} object in the package \code{svfilters.<build>}.
#' @param params  a list of parameters for the amplicon analysis
#' @param transcripts a \code{GRanges} object of the transcripts.  
sv_amplicon_exp <- function(dirs, bviews, grl, amplicon_filters,
                            params=ampliconParams(),
                            transcripts){
  ag_files <- file.path(dirs["2amplicons"],
                        paste0(colnames(bviews), ".rds"))
  result_list <- setNames(vector("list",
                                 length(ag_files)),
                          colnames(bviews))
  for(i in seq_along(result_list)){
    if(file.exists(ag_files[i])){
      ag <- readRDS(ag_files[i])
      result_list[[i]] <- ag
      next()
    }
    ag <- sv_amplicons(bviews[, i], segs=grl[[i]],
                       amplicon_filters, params, transcripts)
    result_list[[i]] <- ag
    saveRDS(ag, file=ag_files[i])
  }
  result_list
}



#' @export
recurrentAmplicons <- function(tx, grl, maxgap=5e3){
  ## tx is big.  Make this smaller as a first step
  tx <- subsetByOverlaps(tx, reduce(unlist(grl), min.gapwidth=maxgap))
  gene_list <- split(tx, tx$gene_name)
  ## remove any element that has multiple chromosomes
  nchroms <- sapply(gene_list, function(g) length(unique(chromosome(g))))
  gene_list <- gene_list[nchroms == 1]
  ## A gene can have multiple entries.  Try reduce
  gene_list <- GRangesList(lapply(gene_list, reduce, min.gapwidth=maxgap))
  ##
  ## Add back the gene name
  ##
  el <- elementNROWS(gene_list)
  tx <- unlist(gene_list)
  tx$gene_name <- rep(names(gene_list), el)
  ## ensure that 2 amplicons for a subject hitting a gene are only counted once
  is_overlap_list <- lapply(grl, function(gr, tx, maxgap) {
    overlapsAny(tx, gr, maxgap=maxgap)
  }, maxgap=maxgap, tx=tx)
  is_overlap <- do.call(cbind, is_overlap_list)
  cnts <- rowSums(is_overlap)
  tx <- tx[cnts > 1, ]
  is_overlap <- is_overlap[ cnts > 1, ]
  cnts <- cnts[ cnts > 1 ]

  ids <- apply(is_overlap, 1, function(is_amplicon, id) {
    paste(id[is_amplicon], collapse=",")
  }, id=gsub(".bam", "", colnames(is_overlap)))
  result <- data.frame(gene = tx$gene_name, freq=as.integer(cnts), id=ids)
  ##
  ## Gene coordinates
  ##
  tx2 <- tx[tx$gene_name %in% result$gene]
  tx2.list <- GRangesList(sapply(split(tx2, tx2$gene_name), reduce))
  tx2 <- unlist(tx2.list)
  tx2 <- tx2[result$gene]
  stopifnot(identical(names(tx2), as.character(result$gene)))
  result$chr <- chromosome(tx2)
  result$start <- start(tx2)
  result$end <- end(tx2)
  result <- result[, c("gene", "chr", "start", "end", "freq", "id")]
  rownames(result) <- NULL
  result[order(result$freq, decreasing=TRUE), ]
}

#' Table of recurrent drivers
#'
#' Given a GRangesList of amplicons or deletions (each element of list is a
#' sample), tabulate the frequency of driver deletions or amplifications.
#' 
#' @param grl a \code{GRangesList} organized by sample
#' @param transcripts a \code{GRanges} object of transcripts as provided by the
#'   \code{svfilters.<ucsc_build>} packages
#' @param split A character string indicating the string used to separate
#'   drivers. By default, we assume the drivers associated with a given deletion
#'   or amplicon are separated by a comma followed by a space.
#' @return a \code{data.frame} of recurrent drivers
#' @export
recurrentDrivers <- function(grl, transcripts, split=", "){
  driver_list <- sapply(grl, function(g){
    dr <- as.character(g$driver)
    if(length(dr) == 0) return(NULL)
    dr <- dr[!is.na(dr)]
    dr <- unlist(strsplit(dr, split))
    unique(dr)
  })
  driver_list2 <- driver_list[ elementNROWS(driver_list) > 0 ]
  df <- data.frame(id=rep(names(driver_list2), elementNROWS(driver_list2)),
                   gene=unlist(driver_list2))
  dflist <- split(df, df$gene)
  ids <- sapply(sapply(dflist, "[[", "id"), paste, collapse=",")
  df <- data.frame(gene=names(dflist),
                   id=ids,
                   freq=elementNROWS(dflist))
  tx <- transcripts[transcripts$gene_name %in% df$gene]
  txlist <- split(tx, tx$gene_name)
  tx <- unlist(reduce(txlist, min.gapwidth=2e3))
  df$chr <- chromosome(tx)
  df$start <- start(tx)
  df$end <- end(tx)
  df <- df[order(df$freq, decreasing=TRUE), ]
  df
}
