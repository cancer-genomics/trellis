#' @include help.R
NULL

setGeneric("clusterReadPairs", function(object) standardGeneric("clusterReadPairs"))
setGeneric("reviseJunction", function(object) standardGeneric("reviseJunction"))
setGeneric("variant<-", function(object, value) standardGeneric("variant<-"))

setClass("LeftAlignmentPairs", contains="GAlignmentPairs")

setValidity("LeftAlignmentPairs", function(object){
  msg <- TRUE
  if(!all(firstIsLeft(object))) return("first read in pair must be the left-most read")
##  index <- order(start(first(object)))
##  if(!identical(index, seq_along(object))) return("left reads must be ordered by physical position")
  msg
})

setReplaceMethod("variant", c("StructuralVariant", "GRanges"),
                 function(object, value){
                   object@variant <- value
                   object
                 })

setMethod("clusterReadPairs", "LeftAlignmentPairs", function(object){
  boundary_starts <- start(first(object))
  boundary_ends <- start(last(object))
  cumsum(c(0, diff(boundary_starts) > 200 | diff(boundary_ends) > 200))
})

setMethod("combine", c("StructuralVariant", "StructuralVariant"),
          function(x, y, ...){
            cnv <- c(variant(x), variant(y))
            prp <- c(proper(x), proper(y))
            irp <- c(improper(x), improper(y))
            cn <- c(copynumber(x), copynumber(y))
            calls <- c(calls(x), calls(y))

            prp_index <- initializeProperIndex2(cnv, prp, zoom.out=1)
            irp_index1 <- initializeImproperIndex2(cnv, irp, DeletionParam())
            irp_index2 <- updateImproperIndex2(cnv, irp, maxgap=2e3)
            irp_index3 <- .match_index_variant(irp_index1, cnv, irp_index2)
            sv <- StructuralVariant(variant=cnv,
                                    proper=prp,
                                    improper=irp,
                                    copynumber=cn,
                                    calls=calls,
                                    index_proper=prp_index,
                                    index_improper=irp_index3)
            sv
          })

setMethod("reviseJunction", "StructuralVariant", function(object){
  v <- variant(object)
  if(length(v) > 1) stop("method only defined for length-one StructuralVariant")
  irp <- improper(object)
  if(length(irp) < 5) return(variant(object))
  irp2 <- leftAlwaysFirst(irp)
  left_border <- end(first(irp2))
  right_border <- start(last(irp2))
  d <- right_border-left_border
  if(any(d > 500)){
    st <- max(left_border[d > 500])+50
    en <- min(right_border[d > 500])-50
    if(st > en) {
      st <- max(left_border)
      en <- min(right_border)
      if(st > en){
        ##stop("Interval for SV is too small after revision -- start
        ##location is larger than end location")
        ##
        ## Keep original boundary
        ##
        st <- start(v)
        en <- end(v)
      }
    }
  } else return(v)
  v <- GRanges(seqnames(v), IRanges(st, en))
  v
})


#' Assess whether a region can be attributable to a germline CNV or
#' artifact seen in germline samples
#'
#' In our application, \code{all_filters} is a \code{GRanges} object
#' comprised of germline CNVs and sequence filters.  The latter is a
#' collection of 1kb bins (reduced) that have low mappability (< 0.75)
#' and/or low GC content (< 0.1).  For each interval in the candidate
#' somatic \code{GRanges} object (\code{g}), we compute the fraction
#' of the interval spanned by the germline filters.  If the fraction
#' of the candidate somatic variant spanned by the germline filter is
#' less than \code{max_proportion_in_filter} AND the width not spanned
#' by the germline filters is less than \code{min_width}, this
#' function evaluates to \code{TRUE}.
#'
#' @return a locical vector having the same length as \code{g}
#' 
#' @seealso See \code{DeletionParam} for default values of
#'   \code{max_proportion_in_filter} and \code{min_width}.
#' 
#' @param g a \code{GRanges} object (e.g., candidate deletions)
#' @param all_filters a \code{GRanges} object (e.g., deletions,
#'   amplifications, and outliers identified in the germline)
#' @param param an instance of \code{DeletionParam}
#' @export
isNotGermline <- function(g, all_filters, param=DeletionParam()){
  p <- intOverWidth(g, all_filters)
  w_not_spanned <- width(g) - p*width(g)
  (p < param@max_proportion_in_filter) &
    (w_not_spanned > param@min_width)
}

isLargeHemizygous <- function(g, param=DeletionParam()){
  (g$seg.mean > homozygousThr(param)) &
    (width(g) > maximumWidth(param))
}

granges_copynumber2 <- function(gr, bins){
  hits <- findOverlaps(gr, bins)
  if (!any(duplicated(names(gr))) && !is.null(names(gr))) {
    split_by <- factor(names(gr)[queryHits(hits)], levels = names(gr))
  }
  else {
    split_by <- queryHits(hits)
  }
  fc_context <- tapply(bins$log_ratio[subjectHits(hits)], queryHits(hits), median)
  as.numeric(fc_context)
}

SVFilters <- function(sv, all_filters, bins, zoom.out=1, param){
  evaluate_context <- evaluateContext(variant(sv), bins)
  if(evaluate_context){
    svcontext <- expandGRanges(variant(sv),
                               10*zoom.out*width(variant(sv))) ## 5percent window
    fc_context <- granges_copynumber2(svcontext, bins)
    ##fc <- 2^(fc_context-copynumber(sv))
    fc <- 2^(copynumber(sv) - fc_context)
  }
  broad_hemizygous <- copynumber(sv) > homozygousThr(param) &
    width(variant(sv)) > maximumWidth(param)
  frac <- intOverWidth(variant(sv), all_filters)
  w <- widthNotSpannedByFilter(variant(sv), all_filters)
  is_dup <- duplicated(variant(sv))
  if(evaluate_context){
    sv <- sv[frac < 0.75 & w > 2e3 & fc < 0.7 & !is_dup & !broad_hemizygous]
  } else{
    sv <- sv[frac < 0.75 & w > 2e3 & !is_dup & !broad_hemizygous]
  }
  sv
}

#' Identify focal, somatic hemizygous deletions and somatic homozygous
#' deletions
#'
#' Segmentation of bin-level normalized log2 ratios yields a set of genomic
#' intervals each having an associated mean. We define candidate somatic
#' deletions as those having a segment mean less than a specified cutoff (e.g.,
#' the theoretical segment mean of a hemizygous deletion on log2 scale). For
#' each candidate CNV, we assess whether the variant can be explained by events
#' identified in germline samples processed in the same batch, or by sequence
#' artifacts (low mappability and/or GC). To help ensure that the deletion is
#' focal (not a small deletion embedded in a larger deletion), we (i) verify
#' that the candidate CNV (\code{cnv}) has a mean log ratio less than the mean
#' log ratio of the neighboring segments by an amount of approx. -0.5
#' [log2(0.7)] and (ii) assess whether the interval could plausibly be a large
#' hemizygous deletion. Candidate somatic deletions identified by this function
#' have the following properties: (i) not germline, (ii) unlikely to be large
#' hemizygous deletions, (iii) the difference in mean of the segment and
#' neighboring segments is less than -0.5, and (iv) have a width of at least
#' 2kb.
#'
#' @return a named \code{GRanges} object.  The names are given by
#'   \code{paste0("sv", seq_along(g))} where g is the \code{GRanges}
#'   object.
#'
#' @param preprocess a list of preprocessed data. See \code{preprocessData}
#' @param germline_filters A \code{GRanges} object (e.g., germline
#'   CNVs and regions of low sequence quality)
#' @param param A \code{DeletionParam} object
#' @seealso \code{\link{preprocessData}}
#' @export
germlineFilters <- function(preprocess, germline_filters, param=DeletionParam()){
  cnv <- preprocess$segments
  bins <- preprocess$bins
  if(missing(germline_filters)){
    germline_filters <- genomeFilters(preprocess$genome)
  }
  if(!is.null(germline_filters)){
    not_germline <- isNotGermline(cnv, germline_filters, param)
  } else {
    ## assume for now that none of the cnvs are germline
    not_germline <- rep(TRUE, length(cnv))
  }
  ##
  ## For unit testing, we may not be able to evaluate the context
  ##
  evaluate_context <- evaluateContext(cnv, bins)
  is_big <- isLargeHemizygous(cnv, param)
  if(evaluate_context){
    egr <- expandGRanges(cnv, 15*width(cnv))
    fc_context <- granges_copynumber2(egr, bins)
    fc <- 2^(cnv$seg.mean-fc_context)
    select <- !is_big & not_germline & fc < 0.7
  }  else{
    select <- !is_big & not_germline
  }
  cnv <- cnv[select]
  if(length(cnv) == 0) {
    return(cnv)
  }
  cnvr <- reduce(cnv)
  K <- width(cnvr) > 2e3
  cnvr <- cnvr[K]
  cnv <- cnv[K]
  names(cnv) <- paste0("sv", seq_along(cnv))
  cnv
}

evaluateContext <- function(g, bins){
  width.cnv <- sum(as.numeric(width(g)))
  width.bins <- sum(as.numeric(width(bins)))
  width.bins/width.cnv > 15
}

initializeImproperIndex <- function(sv, param){
  irp <- sv@improper
  hits <- findOverlaps(variant(sv), irp, minimumGapWidth(param))
  subj_hits <- subjectHits(hits)
  index_improper <- split(subj_hits, names(variant(sv))[queryHits(hits)])
  index_improper
}

initializeImproperIndex2 <- function(gr, improper_rp, param=DeletionParam()){
  hits <- findOverlaps(gr, improper_rp, minimumGapWidth(param))
  subj_hits <- subjectHits(hits)
  index_improper <- split(subj_hits, names(gr)[queryHits(hits)])
  result <- setNames(vector("list", length(gr)), names(gr))
  result[names(index_improper)] <- index_improper
  result
}

initializeImproperIndex3 <- function(sv, param){
  initializeImproperIndex2(variant(sv), sv@improper, param)
}

updateImproperIndex <- function(sv, maxgap=2e3){
  left_boundary <- resize(variant(sv), width=2)
  right_boundary <- resize(variant(sv), width=2, fix="end")
  irp <- sv@improper
  hitsLeft <- findOverlaps(left_boundary, irp, maxgap=maxgap)
  hitsRight <- findOverlaps(right_boundary, irp, maxgap=maxgap)

  index_left <- split(subjectHits(hitsLeft),
                      names(left_boundary)[queryHits(hitsLeft)])
  index_right <- split(subjectHits(hitsRight),
                       names(right_boundary)[queryHits(hitsRight)])

  ## concatenate indices for each element of the index list
  index_all <- setNames(vector("list", length(sv)), names(variant(sv)))
  index_right2 <- index_left2 <- index_all
  i <- match(names(index_left), names(index_all))
  index_left2[i] <- index_left
  i <- match(names(index_right), names(index_all))
  index_right2[i] <- index_right
  updated_index_list <- mapply(function(x,y) unique(intersect(x,y)),
                               index_left2, index_right2)

  ## some of the elements in this list may be NULL
  ##
  ## Assess whether there are any rearranged read pair clusters
  ## 1. order read pairs by left most alignment
  irp_id <- lapply(updated_index_list, function(i, irp) names(irp)[i], irp=irp)
  lrp <- leftAlwaysFirst(irp)##, names(irp)
  lrplist <- lapply(irp_id, function(id, lrp) lrp[names(lrp) %in% id], lrp=lrp)
  ## 2. cluster the read pairs for each element
  cl_list <- lapply(lrplist, clusterReadPairs)
  ## 3. identify id of cluster with most read pairs
  clid <- lapply(cl_list, function(cl) names(which.max(table(cl))))
  ## 4. extract the names of the improper read pairs that correspond
  ## to the above cluster
  irp_id2 <- mapply(function(lrp, clid, cl) names(lrp)[cl==clid],
                    lrp=lrplist, clid=clid, cl=cl_list)
  ## Update the index a second time to include only those improper
  ## read pairs belonging to the cluster
  index_irp <- lapply(irp_id2, function(id, irp) match(id, names(irp)), irp=irp)
  ## the NULLs are now converted to NAs.
  ## Subsetting a vector by NULL returns an empty vector. Convert NAs back to nulls
  na_index <- which(sapply(index_irp, function(x) any(is.na(x))))
  for(j in na_index){
    index_irp[j] <- list(NULL)
  }
  index_irp
}

.hits <- function(query, subject, ...){
  hits <- findOverlaps(query, subject, ...)
  setNames(subjectHits(hits), names(query)[queryHits(hits)])
}

.hits_union <- function(i1, i2){
  nms <- unique(c(names(i1), names(i2)))
  x <- lapply(nms, function(id, i1, i2){
    c(i1[names(i1)==id], i2[names(i2)==id])
  }, i1=i1, i2=i2)
  names(x) <- nms
  x
}

.hits_intersection <- function(i1, i2){
  nms <- unique(c(names(i1), names(i2)))
  x <- lapply(nms, function(id, i1, i2){
    intersect(i1[names(i1)==id], i2[names(i2)==id])
  }, i1=i1, i2=i2)
  names(x) <- nms
  x
}

## Returns a list of indices
## - list is of the same length as the number of variants (length of gr)
##
## - each element contains a vector of indices into the improper read pairs
##  (irp) object
updateImproperIndex2 <- function(gr, irp, maxgap=2e3){
  left_boundary <- resize(gr, width=2)
  right_boundary <- resize(gr, width=2, fix="end")
  if(FALSE){
    ##irp <- sv@improper
    fst <- granges(first(irp))
    lst <- granges(last(irp))
    i1 <- .hits(left_boundary, fst, maxgap=maxgap)
    i4 <- .hits(right_boundary, lst, maxgap=maxgap)
    hits1 <- .hits_intersection(i1, i4)

    i2 <- .hits(left_boundary, lst, maxgap=maxgap)
    i3 <- .hits(right_boundary, fst, maxgap=maxgap)
    hits2 <- .hits_intersection(i2, i3)
    index_list <- mapply(function(x,y) unique(c(x, y)),
                         hits1, hits2)
    return(index_list)
  }
  hitsLeft <- findOverlaps(left_boundary, irp, maxgap=maxgap)
  hitsRight <- findOverlaps(right_boundary, irp, maxgap=maxgap)

  index_left <- split(subjectHits(hitsLeft),
                      names(left_boundary)[queryHits(hitsLeft)])
  index_right <- split(subjectHits(hitsRight),
                       names(right_boundary)[queryHits(hitsRight)])

  ## concatenate indices for each element of the index list
  index_all <- setNames(vector("list", length(gr)), names(gr))
  index_right2 <- index_left2 <- index_all
  i <- match(names(index_left), names(index_all))
  index_left2[i] <- index_left
  i <- match(names(index_right), names(index_all))
  index_right2[i] <- index_right
  updated_index_list <- mapply(function(x,y) unique(intersect(x,y)),
                               index_left2, index_right2)

  ## some of the elements in this list may be NULL
  ##
  ## Assess whether there are any rearranged read pair clusters
  ## 1. order read pairs by left most alignment
  irp_id <- lapply(updated_index_list, function(i, irp) names(irp)[i], irp=irp)
  lrp <- leftAlwaysFirst(irp) ##, names(irp)
  lrplist <- lapply(irp_id, function(id, lrp) lrp[names(lrp) %in% id], lrp=lrp)
  ## 2. cluster the read pairs for each element
  cl_list <- lapply(lrplist, clusterReadPairs)
  ## 3. identify id of cluster with most read pairs
  clid <- lapply(cl_list, function(cl) names(which.max(table(cl))))
  ## 4. extract the names of the improper read pairs that correspond
  ## to the above cluster
  irp_id2 <- mapply(function(lrp, clid, cl) names(lrp)[cl==clid],
                    lrp=lrplist, clid=clid, cl=cl_list)
  ## Update the index a second time to include only those improper
  ## read pairs belonging to the cluster
  index_irp <- lapply(irp_id2, function(id, irp) match(id, names(irp)), irp=irp)
  ## the NULLs are now converted to NAs.
  ## Subsetting a vector by NULL returns an empty vector. Convert NAs back to nulls
  na_index <- which(sapply(index_irp, function(x) any(is.na(x))))
  for(j in na_index){
    index_irp[j] <- list(NULL)
  }
  index_irp
}

## initializeProperIndex <- function(sv, zoom.out=1){
##   proper_rp <- sv@proper
##   g <- expandGRanges(variant(sv), zoom.out*width(variant(sv)))
##   hits <- findOverlaps(g, proper_rp)
##   index_proper <- split(subjectHits(hits), names(variant(sv))[queryHits(hits)])
##   index_proper
## }

initializeProperIndex2 <- function(gr, proper_rp, zoom.out=1){
  g <- expandGRanges(gr, zoom.out*width(gr))
  hits <- findOverlaps(g, proper_rp)
  index_proper <- split(subjectHits(hits), names(gr)[queryHits(hits)])
  result <- setNames(vector("list", length(gr)), names(gr))
  result[names(index_proper)] <- index_proper
  result
}

initializeProperIndex3 <- function(sv, zoom.out=1){
  initializeProperIndex2(variant(sv), sv@proper, zoom.out=zoom.out)
}

improperRP <- function(gr, irp, param=DeletionParam()){
  irp <- updateObject(irp)
  si <- seqinfo(gr)
  sl <- seqlevels(si)
  seqlevels(irp, pruning.mode="coarse") <- sl
  d <- abs(start(first(irp)) - start(last(irp)))
  irp <- irp[d < maximumWidth(param)]
  seqlevelsStyle(irp) <- seqlevelsStyle(si)
  irp <- irp[chromosome(first(irp)) == chromosome(last(irp))]
  hits <- findOverlaps(gr, irp, maxgap=minimumGapWidth(param))
  irp <- irp[unique(subjectHits(hits))]
  if(length(irp) > 0){
    names(irp) <- paste0("i", seq_along(irp))
  }
  irp  
}

addImproperReadPairs2 <- function(gr, aview, param=DeletionParam()){
  irp <- readRDS(improperPaths(aview))
  irp <- updateObject(irp)
  si <- seqinfo(aview)
  sl <- seqlevels(si)
  irp <- keepSeqlevels(irp, sl, pruning.mode="coarse")
  d <- abs(start(first(irp)) - start(last(irp)))
  irp <- irp[d < maximumWidth(param)]
  seqlevelsStyle(irp) <- seqlevelsStyle(si)
  irp <- irp[chromosome(first(irp)) == chromosome(last(irp))]
  hits <- findOverlaps(gr, irp, minimumGapWidth(param))
  irp <- irp[unique(subjectHits(hits))]
  if(length(irp) > 0){
    names(irp) <- paste0("i", seq_along(irp))
  }
  irp
}

.match_index_variant <- function(index_improper, gr, value){
  i <- match(names(value), names(gr))
  index_improper[i] <- value
  index_improper
}


#' Creates a StructuralVariant object, linking genomic intervals
#' (deletions) to read pairs supporting the rearrangement
#'
#' Identifies somatic, focal hemizygous and homozygous deletions, then
#' links these deletions with properly and improperly paired reads.
#'
#' @seealso See \link{germlineFilters} for identifiying somatic deletions.
#' @param aview An \code{AlignmentViews} object
#' @param pview A \code{PreprocessViews} object
#' @param cnv A \code{GRanges} object of candidate somatic deletions
#' @param gr_filters A \code{GRanges} object of germline CNVs and
#'   sequence-based filters identified by low mappability and/or GC
#'   content
#' @param param A \code{DeletionParam} object
deletion_call <- function(preprocess, 
                          gr_filters,
                          param=DeletionParam()){
  if(missing(gr_filters)){
    gr_filters <- genomeFilters(preprocess$genome)
  }
  cnv <- germlineFilters(preprocess, gr_filters, param)
  if(length(cnv) == 0) {
    return(StructuralVariant())
  }
  thr <- homozygousThr(param)
  cncalls <- ifelse(cnv$seg.mean < thr, "homozygous", "hemizygous")
  ##prp <- properReadPairs(bam_path=preprocess$bam.file,
  ##                       gr=cnv, param=param)
  rps <- preprocess$read_pairs
  prp <- prp[["proper_del"]]
  prp_index <- initializeProperIndex2(cnv, prp, zoom.out=1)
  ## change name to filterImproperReadPairs
  ##irp <- addImproperReadPairs2(cnv, aview, param=param)
  improper_rp <- rps[["improper"]]
  irp <- improperRP(cnv, preprocess$improper_rp, param=param)
  irp_index1 <- initializeImproperIndex2(cnv, irp, param)
  irp_index2 <- updateImproperIndex2(cnv, irp, maxgap=2e3)
  irp_index3 <- .match_index_variant(irp_index1, cnv, irp_index2)
  ##irp <- irp[validPairForDeletion(irp)]
  ##irp_index <- initializeImproperIndex2(cnv, irp, param)
  sv <- StructuralVariant(variant=cnv,
                          proper=prp,
                          improper=irp,
                          copynumber=cnv$seg.mean,
                          calls=cncalls,
                          index_proper=prp_index,
                          index_improper=irp_index3)

}

checkHemizygousSpanningHomozygous <- function(object, hits, param, bins){
  spanning_variant <- variant(object)[subjectHits(hits)]
  cncall <- setNames(calls(object), names(variant(object)))
  ## check the part of the interval not containing the homozygous deletion
  portion_notspanning <- setdiff(spanning_variant, variant(object)[queryHits(hits)])
  means <- granges_copynumber2(portion_notspanning, bins)
  means <- sum(means*width(portion_notspanning))/sum(width(portion_notspanning))
  if(means > homozygousThr(param)){
    lirp <- elementNROWS(indexImproper(object))
    L <- lirp[subjectHits(hits)]
    tmp <- ifelse(L <= 4, "hemizygous", "hemizygous+")
    cncall[names(spanning_variant)] <- tmp
  } else cncall[names(spanning_variant)] <- "homozygous"
  cncall[names(spanning_variant)]
}

isHemizygousDeletion <- function(object, param, bins){
  is_hemizygous <- copynumber(object) >= homozygousThr(param) |
    is.na(copynumber(object))
  names(is_hemizygous) <- names(variant(object))
  ## If a hemizygous region contains a homozygous region, the fold
  ## change will be very small but it cwould still be hemizygous...
  ## called "homozygous |-----     ---|  fold change small because of gap.
  ## homozygous               |    |
  not_hemizygous <- !is_hemizygous
  checkhom <- object[not_hemizygous]
  hits <- findOverlaps(variant(checkhom), type="within")
  hits <- hits[subjectHits(hits) != queryHits(hits)]
  ## check whether spanning segment is really hemizygous
  if(length(hits) > 0){
    for(j in seq_len(length(hits))){
      hitsj <- hits[j]
      cncall <- checkHemizygousSpanningHomozygous(checkhom, hitsj, param, bins)
      is_hemizygous[names(cncall)] <- length(grep("hemizygous", cncall)) > 0
    }
  }
  is_hemizygous
}


#' Counts the number of improper read pairs supporting a deletion
#'
#' If the number of improper read pairs supporting a deletion exceeds
#' \code{minFlankingHemizgyous(param)}, a deletion is designated as
#' "hemizygous+" if hemizygous and "homozygous+" if homozygous.
#'
#' @return a character-vector of deletion calls
#' @export
#' @param sv a \code{StructuralVariant} object
#' @param param a \code{DeletionParam} object
#' @param bins a \code{GRanges} object containing the bins used in preprocessing
#'   with mcols value 'log_ratio'
rpSupportedDeletions <- function(sv, param, bins){
  ## NA means that there were no queryRanges in the view -- i.e., all bins were masked
  is_hemizygous <- isHemizygousDeletion(sv, param, bins)
  cncalls <- ifelse(is_hemizygous, "hemizygous", "homozygous")
  number_improper <- elementNROWS(sapply(sv, improper))
  MIN <- param@nflanking_hemizygous
  cncalls[number_improper >= MIN] <- paste0(cncalls[number_improper >= MIN], "+")
  as.character(cncalls)
}


.safelyChangeStart <- function(gr1, gr2){
  if(start(gr2) < end(gr1)){
    return(start(gr2))
  }
  return(start(gr1))
}

.safelyChangeEnd <- function(gr1, gr2){
  if(end(gr2) > start(gr1)){
    return(end(gr2))
  }
  end(gr1)
}

reviseDeletionBorders <- function(object){
  gr <- variant(object)
  for(i in seq_along(gr)){
    X <- object[i]
    g <- reviseJunction(X)
    start(gr)[i] <- .safelyChangeStart(gr[i], g)
    end(gr)[i] <- .safelyChangeEnd(gr[i], g)
  }
  gr
}


.homozygousBorders <- function(object, svs){
  v <- variant(object)
  lrp <- leftAlwaysFirst(proper(object))
  d <- c(0, cumsum(diff(start(first(lrp))) > 500))
  ##
  ## if there are only two clusters, revise the boundaries if close
  ## (if the current boundary is slightly too big, the
  ## findOverlaps(, within) below will not have any hits)
  ##
  if(length(table(d)) == 2) {
    st <- max(end(last(lrp))[d==0])+1
    en <- min(start(first(lrp))[d==1])-1
    ## revise boundaries if 'close' to original
    delta1 <- abs(st-start(v)) < 1000
    delta2 <- abs(en-end(v)) < 1000
    if(delta1) start(v) <- st
    if(delta2) end(v) <- en
  }
  ## now, check whether the variant is within the gap
  segs <- readPairsAsSegments(proper(object))
  g <- gaps(segs)
  hits <- findOverlaps(v, g, type="within")
  if(length(hits) > 0){
    d <- g[subjectHits(hits)]
    start(svs) <- .safelyChangeStart(svs, d)
    end(svs) <- .safelyChangeEnd(svs, d)
  } else svs <- v
  svs
}

##
## resolving homozygous boundaries
##  truth       :  -----|      |---------
##  segmentation:  -------|   |----------
homozygousBorders <- function(object, svs){
  index <- which(calls(object)=="homozygous")
  for(i in index){
    svs[i] <- .homozygousBorders(object[i], svs[i])
  }
  svs
}

addVariant <- function(v, object, cn, cncall, param){
  svtmp <- StructuralVariant(variant=v,
                             improper=object@improper,
                             proper=object@proper,
                             copynumber=cn,
                             calls=cncall)
  indexImproper(svtmp) <- initializeImproperIndex3(svtmp, param)
  indexProper(svtmp) <- initializeProperIndex3(svtmp, zoom.out=1)
  indexImproper(svtmp) <- updateImproperIndex(svtmp, maxgap=1e3)

  cnvs <- c(v, granges(variant(object)))
  ids <- paste0("sv", seq_along(cnvs))
  cnvs <- setNames(cnvs, ids)
  cncalls <- as.character(c(cncall, calls(object)))
  cns <- as.numeric(c(cn, copynumber(object)))
  index_proper <- setNames(c(indexProper(svtmp), indexProper(object)), ids)
  index_improper <- setNames(c(indexImproper(svtmp), indexImproper(object)), ids)
  sv2 <- StructuralVariant(cnvs,
                           calls=cncalls,
                           copynumber=cns,
                           proper=object@proper,
                           improper=object@improper,
                           index_proper=index_proper,
                           index_improper=index_improper)
  sv2 <- sv2[!duplicated(cnvs)]
  sv2
}

allGapsInReadPairs <- function(proper.readpairs){
  segs <- readPairsAsSegments(proper.readpairs)
  g <- gaps(segs)
  g <- g[width(g) > 2e3]
  g
}

gapsInHemiDel <- function(proper.readpairs, deletion.gr){
  all.gaps <- allGapsInReadPairs(proper.readpairs)
  expanded.deletion <- expandGRanges(deletion.gr, 1e3)
  gaps.in.hemi.del <- subsetByOverlaps(all.gaps, expanded.deletion, type="within")
  if(length(gaps.in.hemi.del) > 0){
    si <- seqinfo(deletion.gr)
    seqlevels(gaps.in.hemi.del, pruning.mode="coarse") <- seqlevels(si)
    seqinfo(gaps.in.hemi.del) <- si    
  }
  gaps.in.hemi.del
}

spansHomozygousDel <- function(sv){
  hits <- homozygousHits(sv)
  length(hits > 0)
}

addVariant2 <- function(v, object, cn, cncall, param){
  v$seg.mean <- cn
  v$sample <- variant(object)$sample[1]
  svtmp <- StructuralVariant(variant=v,
                             improper=object@improper,
                             proper=object@proper,
                             copynumber=cn,
                             calls=cncall)
  indexImproper(svtmp) <- initializeImproperIndex3(svtmp, param)
  indexProper(svtmp) <- initializeProperIndex3(svtmp, zoom.out=1)
  indexImproper(svtmp) <- updateImproperIndex(svtmp, maxgap=1e3)

  ##cnvs <- c(v, granges(variant(object)))
  cnvs <- c(v, variant(object))
  ids <- paste0("sv", seq_along(cnvs))
  cnvs <- setNames(cnvs, ids)
  cncalls <- as.character(c(cncall, calls(object)))
  cns <- as.numeric(c(cn, copynumber(object)))
  index_proper <- setNames(c(indexProper(svtmp), indexProper(object)), ids)
  index_improper <- setNames(c(indexImproper(svtmp), indexImproper(object)), ids)
  sv2 <- StructuralVariant(cnvs,
                           calls=cncalls,
                           copynumber=cns,
                           proper=object@proper,
                           improper=object@improper,
                           index_proper=index_proper,
                           index_improper=index_improper)
  sv2 <- sv2[!duplicated(cnvs)]
  sv2
}


## Check for overlapping hemizygous deletions that create a homozygous deletion
##  truth1       :  -----|      |--------------
##  truth2       :|---------|        |---------
##  homozyg gap  :         ->   <-
.hemizygousBorders <- function(del.gr, sv, param){
  gaps.in.hemi.del <- gapsInHemiDel(proper(sv), del.gr)
  if(length(gaps.in.hemi.del) == 0) return(sv)
  ##
  ## the gap is potentially a homozygous deletion formed by overlapping
  ## hemizygous deletions
  ##
  epsilon <- rep(log2(1/50), length(gaps.in.hemi.del))
  sv2 <- addVariant2(gaps.in.hemi.del,
                     object=sv,
                     cn=epsilon,
                     cncall=rep("homozygous", length(gaps.in.hemi.del)),
                     param=param)
  sv2
}

hemizygousBorders <- function(object, param){
  index <- grep("hemizygous", calls(object))
  sv <- object
  hemizygous.gr <- variant(object)[index]
  for(i in seq_along(hemizygous.gr)){
    sv <- .hemizygousBorders(del.gr=hemizygous.gr[i],
                              sv=sv, param=param)
  }
  sv
}


improperReadPairs <- function(aview, gr, param=DeletionParam()){
  irp <- readRDS(improperPaths(aview))
  irp <- updateObject(irp)
  si <- seqinfo(aview)
  sl <- seqlevels(si)
  irp <- keepSeqlevels(irp, sl, pruning.mode="coarse")
  d <- abs(start(first(irp)) - start(last(irp)))
  irp <- irp[d < maximumWidth(param)]
  seqlevelsStyle(irp) <- seqlevelsStyle(si)
  irp <- irp[chromosome(first(irp)) == chromosome(last(irp))]
  hits <- findOverlaps(gr, irp, minimumGapWidth(param))
  irp <- irp[unique(subjectHits(hits))]
  if(length(irp) > 0){
    names(irp) <- paste0("i", seq_along(irp))
  }
  irp
}

.reviseEachJunction <- function(sv, bins, improper_rp, param=DeletionParam()){
  svs <- reviseDeletionBorders(sv)
  ##
  ## for the duplicated ranges, revert back to the original
  ##
  svs[duplicated(svs)] <- variant(sv)[duplicated(svs)]
  variant(sv) <- svs
  copynumber(sv) <- granges_copynumber2(svs, bins)
  variant(sv) <- homozygousBorders(sv, svs)
  is.dup <- duplicated(variant(sv))
  if(any(is.dup)){
    sv <- sv[!is.dup]
  }
  sv <- hemizygousBorders(sv, param)
  irp <- improperRP(variant(sv), improper_rp, param=param)
  improper(sv) <- irp
  indexImproper(sv) <- updateImproperIndex2(variant(sv), irp, maxgap=500)
  calls(sv) <- rpSupportedDeletions(sv, param, bins)
  sort(sv)
}

removeDuplicateIntervals <- function(g, object){
  is.duplicated <- duplicated(g)
  dups <- g[is.duplicated]
  k <- match(names(dups), names(variant(object)))
  object <- object[-k]
  object
}

adjudicateHemizygousOverlap2 <- function(object){
  cncalls <- gsub("\\+", "", calls(object))
  hemi1_or_hemi2 <- cncalls=="hemizygous" | cncalls=="OverlappingHemizygous"
  if(!any(hemi1_or_hemi2)) return(object)
  g <- variant(object)
  if(length(g) < 2) return(object)
  is.duplicated <- duplicated(g)
  if(any(is.duplicated)){
    object <- removeDuplicateIntervals(g, object)
  }
  g <- g[names(g) %in% names(variant(object))]
  if(length(g) < 2) return(object)
  ##
  ## If one interval is contained in another, drop the interval with the least
  ## support.
  ##
  self_hits <- selfHits(g, type="within")
  if(length(self_hits) == 0){
    ## one interval is not contained within another.  Keep both
    ## -----------          -----------
    ## ------           ---------------
    return(object)
  }
  ## Keep the interval with most support
  ## -------------        --------------
  ## ---------------     ---------------
  ##
  k <- match(names(g), names(variant(object)))
  indices <- indexImproper(object[k])
  drop <- which.min(elementNROWS(indices))
  if(length(drop) > 0){
    dropid <- names(g)[drop]
    dropindex <- match(dropid, names(variant(object)))
    object <- object[-dropindex]
  }
  object
}

adjudicateHemizygousOverlap <- function(g, object){
  if(length(g) < 2) return(object)
  is.duplicated <- duplicated(g)
  if(any(is.duplicated)){
    object <- removeDuplicateIntervals(g, object)
  }
  g <- g[names(g) %in% names(variant(object))]
  if(length(g) < 2) return(object)
  ##
  ## If one interval is contained in another, drop the interval with the least
  ## support.
  ##
  hits <- findOverlaps(g, type="within")
  hits <- hits[subjectHits(hits) != queryHits(hits)]
  if(length(hits) == 0){
    ## one interval is not contained within another.  Keep both
    ## -----------          -----------
    ## ------           ---------------
    return(object)
  }
  ## Keep the interval with most support
  ## -------------        --------------
  ## ---------------     ---------------
  ##
  k <- match(names(g), names(variant(object)))
  indices <- indexImproper(object[k])
  drop <- which.min(elementNROWS(indices))
  if(length(drop) > 0){
    dropid <- names(g)[drop]
    dropindex <- match(dropid, names(variant(object)))
    object <- object[-dropindex]
  }
  object
}

removeDuplicateSv <- function(object){
  g <- variant(object)
  is.duplicated <- duplicated(g)
  if(any(is.duplicated)){
    object <- removeDuplicateIntervals(g, object)
  }
  object
}

selfHits <- function(g, ...){
  hits <- findOverlaps(g, ...)
  hits <- hits[subjectHits(hits) != queryHits(hits)]
  hits
}

.merge_partial_overlap <- function(sv){
  gr <- variant(sv)
  cn <- copynumber(sv)
  cn2 <- sum(cn*width(variant(sv)))/sum(width(variant(sv)))
  ##  nms <- names(g)[m]
  ##  names(rr) <- nms[1]
  ##  dropid <- nms[-1]
  ##  dropindex <- match(dropid, names(variant(object)))
  ##  object <- object[-dropindex]
  sv2 <- sv[1]  ## only keep the first
  rgr <- setNames(reduce(gr), names(gr)[1])
  rgr$seg.mean <- cn2
  rgr$sample <- gr$sample[1]
  variant(sv2) <- rgr
  copynumber(sv2) <- cn2
  indexImproper(sv2) <- updateImproperIndex(sv2, maxgap=1e3)
  indexProper(sv2) <- initializeProperIndex3(sv2, zoom.out=1)
  ##  ## update the improperPair index
  ##  v <- variant(object)
  ##  i <- match(names(rr), names(v))
  ##  v[i] <- rr
  ##  variant(object) <- v
  ##  obj <- object[i]
  ##  indexImproper(obj) <- updateImproperIndex(obj, maxgap=1e3)
  ##  indexProper(obj) <- initializeProperIndex3(obj, zoom.out=1)
  ##  copynumber(obj) <- cn2
  calls(sv2) <- "homozygous"
  ##indices <- indexImproper(sv2)
  ##K <- match(names(rr), names(indices))
  ##indices[[K]] <- indexImproper(obj)[[1]]
  ##object@index_improper <- indices
  ##indices <- indexProper(object)
  ##indices[[K]] <- indexProper(obj)[[1]]
  ##object@index_proper <- indices
  sv2
}

## Check whether any homozygous deletions overlap. If overlap, keep the reduced
## interval
adjudicateHomozygousOverlap2 <- function(object){
  gr <- variant(object)
  gr <- gr[calls(object)=="homozygous"]
  if(length(gr) < 2) return(object)
  sv2 <- removeDuplicateSv(object)
  gr <- variant(sv2)
  if(length(gr) < 2) next()
  self_hits <- selfHits(variant(sv2))
  if(length(self_hits) == 0) return(sv2)
  ##
  ## (1)---------   ----------
  ## (2)-------        -------
  ## Remove (1)
  self_hits <- selfHits(variant(sv2), type="within")
  if(length(self_hits) > 0){
    sv2 <- sv2[-queryHits(self_hits)]
  }
  ##
  ## Overlapping
  ## -------------        --------------
  ## ---------------     ---------------
  ##
  ## Take the union.
  hits <- findOverlaps(variant(sv2), reduce(variant(sv2)))
  indices <- split(queryHits(hits), subjectHits(hits))
  indices <- indices[elementNROWS(indices) > 1]
  svnew <- sv2
  for(k in seq_along(indices)){
    m <- indices[[k]]
    sv3 <- .merge_partial_overlap(sv2[m])
    svnew <- combine(svnew[-m], sv3)
  }
  is.dup <- duplicated(variant(svnew))
  svnew <- svnew[!is.dup]
  svnew
}

adjudicateHomozygousOverlap <- function(g, object){
  if(length(g) < 2) return(object)
  is.duplicated <- duplicated(g)
  if(any(is.duplicated))
    object <- removeDuplicateIntervals(g, object)
  g <- g[names(g) %in% names(variant(object))]
  if(length(g) < 2) return(object)
  ##
  ## Can not have overlapping homozygous deletions.
  ##
  hits <- findOverlaps(g)
  hits <- hits[subjectHits(hits) != queryHits(hits)]
  if(length(hits) == 0){
    ## Intervals not overlapping. Keep both
    ## -----------          -----------
    ## ------    ---------------
    return(object)
  }
  ##
  ## (1)---------   ----------
  ## (2)-------        -------
  ## Remove (1)
  hits <- findOverlaps(g, type="within")
  hits <- hits[subjectHits(hits) != queryHits(hits)]
  if(length(hits) > 0){
    dropid <- names(g)[queryHits(hits)]
    object <- object[-match(dropid, names(variant(object)))]
    g <- g[-match(dropid, names(g))]
  }
  ##
  ## Overlapping
  ## -------------        --------------
  ## ---------------     ---------------
  ##
  ## Take the union.
  r <- reduce(g)
  if(identical(length(r), length(g)))    return(object)
  hits <- findOverlaps(g, r)
  indices <- split(queryHits(hits), subjectHits(hits))
  indices <- indices[elementNROWS(indices) > 1]
  for(k in seq_along(indices)){
    m <- indices[[k]]
    rr <- r[k]
    J <- match(names(g)[m], names(variant(object)))
    cn <- copynumber(object)[J]
    cn2 <- sum(cn*width(variant(object))[J])/sum(width(variant(object))[J])

    nms <- names(g)[m]
    names(rr) <- nms[1]
    dropid <- nms[-1]
    dropindex <- match(dropid, names(variant(object)))
    object <- object[-dropindex]
    ## update the improperPair index
    v <- variant(object)
    i <- match(names(rr), names(v))
    v[i] <- rr
    variant(object) <- v
    obj <- object[i]
    indexImproper(obj) <- updateImproperIndex(obj, maxgap=1e3)
    indexProper(obj) <- initializeProperIndex3(obj, zoom.out=1)
    copynumber(obj) <- cn2
    calls(obj) <- "homozygous"

    indices <- indexImproper(object)
    K <- match(names(rr), names(indices))
    indices[[K]] <- indexImproper(obj)[[1]]
    object@index_improper <- indices

    indices <- indexProper(object)
    indices[[K]] <- indexProper(obj)[[1]]
    object@index_proper <- indices
  }
  is.dup <- duplicated(variant(object))
  object <- object[!is.dup]
  object
}

removeSameStateOverlapping2 <- function(sv){
  cncalls <- gsub("\\+", "", calls(sv))
  is.homdel <- cncalls == "homozygous"
  sv.homdel <- sv[is.homdel]
  sv.hemdel <- sv[!is.homdel]
  sv.hemdel <- adjudicateHemizygousOverlap2(sv.hemdel)
  sv.homdel <- adjudicateHomozygousOverlap2(sv.homdel)
  sv <- combine(sv.hemdel, sv.homdel)
  sv <- sort(sv)
  sv
}

removeSameStateOverlapping <- function(sv){
  v <- variant(sv)
  r <- reduce(v)
  if(length(r) == length(v)) return(sv)
  hits <- findOverlaps(v, r)
  hitlist <- split(queryHits(hits), subjectHits(hits))
  hitlist <- hitlist[elementNROWS(hitlist) > 1]
  if(length(hitlist) == 0) return(sv)
  sv2 <- sv
  for(j in seq_along(hitlist)){
    k <- hitlist[[j]]
    gr <- v[k]
    cncalls <- gsub("\\+", "", calls(sv)[k])
    hemi1_or_hemi2 <- cncalls=="hemizygous" | cncalls=="OverlappingHemizygous"
    sv2 <- adjudicateHemizygousOverlap(gr[hemi1_or_hemi2], sv2)
    gr <- gr[cncalls=="homozygous"]
    if(length(gr) < 2) next()
    sv2 <- removeDuplicateSv(gr, sv2)
    gr <- gr[names(gr) %in% names(variant(sv2))]
    if(length(gr) < 2) next()
    if(!anyOverlaps(gr, variant(sv2))) next()
    sv2 <- adjudicateHomozygousOverlap2(tmp, sv2)
  }
  sv2
}

reviseEachJunction <- function(sv, bins, improper_rp, param=DeletionParam()){
  message("Revising junctions...")
  sv <- .reviseEachJunction(sv, bins, improper_rp, param)
  indexProper(sv) <- initializeProperIndex3(sv, zoom.out=1)
  copynumber(sv) <- granges_copynumber2(variant(sv), bins)
  calls(sv) <- rpSupportedDeletions(sv, param, bins)
  sv
}

findSpanningHemizygousDeletion <- function(hits, homdel, irp, object, bins, param){
  K <- queryHits(hits)[1]
  ##homdel <- variant(object)[K]
  if(length(hits) < 5) return(object) ##nothing to do
  j <- subjectHits(hits)
  irpj <- irp[j]
  rpsegs <- readPairsAsSegments(irpj)
  ## do these segments span the putative homozygous deletion?
  hitw <- findOverlaps(homdel, rpsegs, type="within")
  if(length(hitw) < 5) return(object)
  lrp <- leftAlwaysFirst(irpj[subjectHits(hitw)])
  cl <- clusterReadPairs(lrp)
  clid <- names(which.max(table(cl)))
  lrp <- lrp[cl==clid]
  if(length(lrp) < 5) return(object)
  ##
  ## Add hemizygous deletion
  ##
  st <- max(end(first(lrp)))
  en <- min(start(last(lrp)))
  hemdel <- GRanges(seqnames(homdel)[1], IRanges(start=st, end=en))
  hiteq <- findOverlaps(hemdel, homdel, type="equal")
  if(length(hiteq) > 0){
    ##
    ##1. called homdel:  ------------------               ----------------------
    ##2. true   homdel:  -----------------------          ----------------------
    ##
    ##  see if there are gaps in the proper pairs to identify (2)
    ##
    rps <- object@proper
    rps2 <- readPairsAsSegments(rps)
    rps3 <- rps2[overlapsAny(rps2, homdel)]
    g <- gaps(rps3)
    if(length(g)== 0) {
      ## there are no gaps
      stop("there are no gaps")
    }
    g <- g[overlapsAny(g, homdel)]
    g <- g[width(g) >= 2e3]
    g <- GRanges(seqnames(g)[1], IRanges(min(start(g)), max(end(g))))
    seqlevels(g, pruning.mode="coarse") <- seqlevels(homdel)
    seqinfo(g) <- seqinfo(homdel)
    object2 <- addVariant2(v=g,
                          object=object,
                          cn=log2(1/50),
                          cncall="homozygous",
                          param=param)
    ##
    ## The homdel is really hemizygous
    ##
    hemdel <- homdel
    ## Remove the gaps and recompute the foldchange
    portion_notspanning <- setdiff(hemdel, g)
    means <- granges_copynumber2(portion_notspanning, bins)
    means <- sum(means*width(portion_notspanning))/sum(width(portion_notspanning))
    ##
    ## Update the call and the fold change
    ##
    calls(object2)[K] <- "hemizygous+"
    copynumber(object2)[K] <- means
    return(object2)
  }
  seqlevels(hemdel,pruning.mode="coarse") <- seqlevels(homdel)
  seqinfo(hemdel) <- seqinfo(homdel)
  portion_notspanning <- setdiff(hemdel, homdel)
  means <- granges_copynumber2(portion_notspanning, bins)
  means <- sum(means*width(portion_notspanning))/sum(width(portion_notspanning))
  object2 <- addVariant2(v=hemdel,
                         object=object,
                         cn=means,
                         cncall="hemizygous+",
                         param=param)
  indexImproper(object2) <- updateImproperIndex(object2)
  sort(object2)
}

leftHemizygousHomolog <- function(object, bins, param){
  ## #########################################################################
  ##
  ## Identifying overlapping hemizygous deletions
  ##
  ## #########################################################################
  ##
  ## homolog1 ---------|       |-------------
  ## homolog2 -----|        |----------------
  ##
  ## RRP  1:       |--------|
  ## RRP  2:           |------|
  ##
  ## 1.  Find RRP 1 by getting all improper read pairs near the right
  ## end of the homozygous deletion
  ##
  ## order improper RPs by left-most position
  sv_bak <- object
  homindex <- which(calls(object)=="homozygous")
  irp <- leftAlwaysFirst(object@improper)
  ## define rightedge of homozygous deletion
  gr <- variant(object)
  ##homsv <- sv[homindex]
  rightedge <- restrict(gr, start=end(gr)-2e3)
  end(rightedge) <- end(rightedge)+2e3

  hits <- findOverlaps(rightedge, last(irp))
  ## group hits by the homozygous deletion interval
  id <- names(rightedge)[queryHits(hits)]
  hitlist <- split(hits, factor(id,levels=unique(id)))
  hitlist <- hitlist[names(hitlist) %in% names(gr)[homindex]]
  object2 <- object
  if(length(hitlist) > 0){
    ##
    ## Function for single hit object
    ##
    for(k in seq_along(hitlist)){
      ##
      ## object is stable.  object2 is not!
      ##
      homdel <- variant(object)[names(hitlist)[k]]
      object2 <- findSpanningHemizygousDeletion(hits=hitlist[[k]],
                                                homdel=homdel,
                                                irp=irp,
                                                object=object2,
                                                bins=bins,
                                                param=param)
    }
  }
  object2
}

firstIsLeft <- function(galp) start(first(galp)) <= start(last(galp))

leftAlwaysFirst <- function(rp){
  ##is_r1_left <- start(first(rp)) < start(last(rp))
  is_r1_left <- firstIsLeft(rp)
  rp2 <- GAlignmentPairs(first=c(first(rp)[is_r1_left],
                           last(rp)[!is_r1_left]),
                         last=c(last(rp)[is_r1_left],
                           first(rp)[!is_r1_left]),
                         isProperPair=rep(FALSE, length(rp)))
  ids <- c(names(first(rp)[is_r1_left]), names(last(rp)[!is_r1_left]))
  names(rp2) <- ids
  rp2 <- rp2[order(start(first(rp2)))]
  as(rp2, "LeftAlignmentPairs")
}

rightHemizygousHomolog <- function(object, bins, param){
  ## 2.  Find RRP 2 by getting all improper read pairs near the left
  ## end of the homozygous deletion
  homindex <- which(calls(object)=="homozygous")
  gr <- variant(object)
  leftedge <- restrict(gr, end=start(gr)+2e3)
  start(leftedge) <- start(leftedge)-2e3
  irp <- leftAlwaysFirst(object@improper)

  hits <- findOverlaps(leftedge, first(irp))
  ## group hits by the homozygous deletion interval
  id <- names(leftedge)[queryHits(hits)]
  hitlist <- split(hits, factor(id,levels=unique(id)))
  hitlist <- hitlist[names(hitlist) %in% names(gr)[homindex]]
  object2 <- object
  if(length(hitlist) > 0){
    ##
    ## Function for single hit object
    ##
    for(k in seq_along(hitlist)){
      ##
      ## object is stable.  object2 is not!
      ##
      homdel <- variant(object)[names(hitlist)[k]]
      object2 <- findSpanningHemizygousDeletion(hits=hitlist[[k]],
                                                homdel=homdel,
                                                irp=irp,
                                                object=object2,
                                                bins=bins,
                                                param=param)
    }
  }
  object2 <- sort(object2)
  object2
}

refineHomozygousBoundaryByHemizygousPlus <- function(sv){
  hemizygousplus <- variant(sv)[calls(sv)=="hemizygous+"]
  ##homozygous <- variant(sv)[calls(sv)=="homozygous"]
  homozygous <- variant(sv)[grep("homozygous", calls(sv))]
  hits <- findOverlaps(homozygous, hemizygousplus)
  delta1 <- abs(start(homozygous)[queryHits(hits)] - start(hemizygousplus)[subjectHits(hits)]) < 1000
  delta2 <- abs(end(homozygous)[queryHits(hits)] - end(hemizygousplus)[subjectHits(hits)]) < 1000
  if(any(delta1)){
    jj <- subjectHits(hits)[delta1]
    ii <- queryHits(hits)[delta1]
    start(homozygous)[ii] <- start(hemizygousplus)[jj]
  }
  if(any(delta2)){
    jj <- subjectHits(hits)[delta2]
    ii <- queryHits(hits)[delta2]
    end(homozygous)[ii] <- end(hemizygousplus)[jj]
  }
  variant(sv)[grep("homozygous", calls(sv))] <- homozygous
  sv
}

callOverlappingHemizygous <- function(object){
  hits <- findOverlaps(variant(object), type="within")
  hits <- hits[subjectHits(hits) != queryHits(hits)]
  if(length(hits) > 0){
    j <- subjectHits(hits)
    calls(object)[j] <- "OverlappingHemizygous+"
  }
  object
}

removeSameStateOverlapping <- function(sv){
  v <- variant(sv)
  r <- reduce(v)
  if(length(r) == length(v)) return(sv)
  hits <- findOverlaps(variant(sv), r)
  hitlist <- split(queryHits(hits), subjectHits(hits))
  hitlist <- hitlist[elementNROWS(hitlist) > 1]
  if(length(hitlist) == 0) return(sv)
  sv2 <- sv
  for(j in seq_along(hitlist)){
    k <- hitlist[[j]]
    tmp <- v[k]
    cncalls <- gsub("\\+", "", calls(sv)[k])
    hemi1_or_hemi2 <- cncalls=="hemizygous" | cncalls=="OverlappingHemizygous"
    sv2 <- adjudicateHemizygousOverlap(tmp[hemi1_or_hemi2], sv2)
    sv2 <- adjudicateHomozygousOverlap(tmp[cncalls=="homozygous"], sv2)
  }
  sv2
}

widthNotSpannedByFilter <- function(cnv, filters){
  if(length(filters) != length(reduce(filters))) stop("filter GRanges not reduced")
  w <- width(cnv)
  int <- intersect(cnv, filters)
  if(length(int)== 0) return(w)
  int_width <- width(int)
  cnvhits <- findOverlaps(int, cnv)
  width_of_cnv <- width(cnv)[unique(subjectHits(cnvhits))]
  width_of_cnv_intersecting_filters <- sapply(split(int_width[queryHits(cnvhits)],
                                                    subjectHits(cnvhits)), sum)
  width_not_in_filter <- width_of_cnv - width_of_cnv_intersecting_filters
  w[unique(subjectHits(cnvhits))] <- as.numeric(width_not_in_filter)
  w
}



groupSVs <- function(object){
  g <- variant(object)
  g2 <- expandGRanges(g, width(g)*1)
  gr <- reduce(g2)
  hits <- findOverlaps(variant(object), gr)
  j <- subjectHits(hits)
  groupedVariant(object) <- factor(j)
  object
}


allProperReadPairs <- function(sv, param, bfile, zoom.out=1){
  cnv <- variant(sv)
  g <- expandGRanges(cnv, width(cnv)*zoom.out)
  query <- reduce(disjoin(g))
  proper_rp <- readPairsNearVariant(query, bfile)
  sv@proper <- proper_rp
  indexProper(sv) <- initializeProperIndex3(sv, zoom.out=zoom.out)
  sv
}

removeHemizygous <- function(sv){
  is_hemizygous <- calls(sv)=="hemizygous"
  sv <- sv[!is_hemizygous]
  sv
}

revise <- function(sv, bins, param){
  copynumber(sv) <- granges_copynumber2(variant(sv), bins)
  calls(sv) <- rpSupportedDeletions(sv, param=param, bins)
  indexImproper(sv) <- updateImproperIndex(sv, maxgap=500)
  calls(sv) <- rpSupportedDeletions(sv, param, bins)
  sv2 <- leftHemizygousHomolog(sv, bins, param)
  sv3 <- rightHemizygousHomolog(sv2, bins, param)
  calls(sv3) <- rpSupportedDeletions(sv3, param, bins)
  message("Refining homozygous boundaries by spanning hemizygous+")
  sv5 <- refineHomozygousBoundaryByHemizygousPlus(sv3)
  sv6 <- callOverlappingHemizygous(sv5)
  sv7 <- removeSameStateOverlapping2(sv6) 
  sv7
}

#' Creates a StructuralVariant object for a sample
#'
#' Creates a \code{StructuralVariant} object encapsulating information
#' on deletions, including the genomic intervals, proper and improper
#' read pairs in the vicinity, and a classification for the type of
#' deletion based on the preprocessed estimates of read depth and the
#' improper read pair alignments.
#'
#' Note: proper read pairs near a candidate deletion are randomly sampled to
#' reduce the size of the object. For reproducibility, set a seed prior to
#' running this function.
#'
#' @param preprocess a list of preprocessing summaries (see \code{preprocessData})
#' @param gr_filters a \code{GRanges} object of germline filters
#' @param param a \code{DeletionParam} object
#' @seealso \code{\link{preprocessData}}
#' @export
sv_deletions <- function(preprocess,
                         gr_filters,
                         param=DeletionParam()){
  if(missing(gr_filters)){
    gr_filters <- genomeFilters(preprocess$genome)
  }
  segs <- preprocess$segments
  preprocess$segments <- segs[segs$seg.mean < hemizygousThr(param)]
  sv <- deletion_call(preprocess, gr_filters, param)
  ##sv <- deletion_call(bam.file, improper_rp, pview, gr, gr_filters)
  calls(sv) <- rpSupportedDeletions(sv, param=param, bins=preprocess$bins)
  sv <- removeHemizygous(sv)
  improper_rp <- preprocess$read_pairs[["improper"]]
  sv <- reviseEachJunction(sv, preprocess$bins, improper_rp, param)
  sv <- removeHemizygous(sv)
  sv <- revise(sv, bins=preprocess$bins, param=param)
  ## requires bam file
  sv <- finalize_deletions(sv=sv, preprocess,
                           gr_filters=gr_filters,
                           param=param)
  sv
}

setMethod("rename", "StructuralVariant", function(x, ...){
  nms <- paste0("sv", seq_along(variant(x)))
  names(variant(x)) <- nms
  ip <- indexProper(x)
  names(ip) <- nms
  x@index_proper <- ip
  ip <- indexImproper(x)
  names(ip) <- nms
  x@index_improper <- ip
  names(copynumber(x)) <- nms
  x
})

finalize_deletions <- function(sv, preprocess, gr_filters,
                               param=DeletionParam()){
  if(missing(gr_filters)){
    gr_filters <- genomeFilters(preprocess$genome)
  }
  bins <- preprocess$bins
  bam.file <- preprocess$bam.file
  sv <- SVFilters(sv, gr_filters, bins, param=param)
  sv <- groupSVs(sv)
  ## requires bam file
  sv <- allProperReadPairs(sv, param,
                           bfile=bam.file, zoom.out=1)
  if(length(sv@proper) > 25e3){
    proper(sv) <- sv@proper[sample(seq_along(sv@proper), 25e3)]
    indexProper(sv) <- initializeProperIndex3(sv, zoom.out=1)
  }
  sv <- rename(sort(sv))
  sv
}

#' Extract a GRangesList of all identified deletions in a project
#'
#' Deletions are assumed by be saved in a deletion directory by the name of
#' the BAM file.
#'
#' @param path character string providing directory where deletions are saved
#' @param ids a character vector of sample identifiers
#' @return a \code{GRangesList} of deletions
#' @export
listDeletions <- function(path="data/segment/1deletions", ids){
  ##files <- file.path(dp["1deletions"], paste0(ids, ".rds"))
  files <- file.path(path, paste0(ids, ".rds"))
  dels <- lapply(files, readRDS)
  names(dels) <- ids
  g <- lapply(dels, function(x) granges(variant(x)))
  g <- unlist(GRangesList(g))
  g$id <- sapply(strsplit(names(g), "\\."), "[", 1)
  names(g) <- NULL
  g$log_ratio <- unlist(lapply(dels, copynumber))
  g$call <- unlist(lapply(dels, calls))
  grl <- split(g, g$id)
  grl
}

#' @export
reduceTranscripts <- function(tx, grl, maxgap=5e3){
  ## tx is big.  Make this smaller as a first step
  if(!missing(grl)){
    tx <- subsetByOverlaps(tx, reduce(unlist(grl), min.gapwidth=maxgap))
  }
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
  sort(tx)
}

#' Tabulate the frequency a gene has a deletion
#'
#' Multiple hit genes. Deletions do not necessarily have to overlap so long as
#' the two deletions both hit the same gene.
#'
#' 1. list transcripts by gene
#' 2. reduce the deletion granges for each subject to avoid over-counting
#' overlapping hemizygous deletions as 2 hits
#' 3. count overlaps
#' @return a \code{data.frame} of gene names with frequencies
#' @param tx a \code{transcript} object as provided by the \code{svfilters} package
#' @param grl a \code{GRangesList} of deletions
#' 
#' @param maxgap a length-one numeric vector indicating the amount of space
#'   between a deletion and a transcript allowed to consider overlapping. This
#'   argument is passed to \code{overlapsAny}.
#' 
#' @seealso \code{\link[IRanges]{overlapsAny}}
#' @return a \code{GRanges} object. **The coordinates are of the reduced transcript and not the observed deletion.**
#' @export
numberOverlappingDeletions <- function(tx, grl, maxgap=5e3){
  ## ensure that 2 deletions for a subject hitting a gene are only counted once
  is_overlap_list <- lapply(grl, function(gr, tx, maxgap) {
    overlapsAny(tx, gr, maxgap=maxgap)
  }, maxgap=maxgap, tx=tx)
  is_overlap <- do.call(cbind, is_overlap_list)
  frequency <- rowSums(is_overlap)
  frequency
}

#' @export
numberSpanningAmplicons <- function(tx, grl, maxgap=5e3){
  ## ensure that 2 deletions for a subject hitting a gene are only counted once
  is_overlap_list <- lapply(grl, function(gr, tx, maxgap) {
    overlapsAny(tx, gr, type="within")
  }, maxgap=maxgap, tx=tx)
  is_overlap <- do.call(cbind, is_overlap_list)
  frequency <- rowSums(is_overlap)
  frequency
}




#' Extract read pairs from a \code{StructuralVariant} object
#'
#' This function extracts all read pairs from a StructuralVariant object and
#' then performs a filtering step. The filtering step thins the number of read
#' pairs with insert size less than 3kb.
#'
#' @param object  a \code{StructuralVariant} object
#' @param max.out length-one integer vector: the maximum number of read pairs to return
#' @export
thinReadPairs <- function(object, max.out=25e3){
  rp <- readPairs(object)
  L <- length(rp)
  if(L > max.out){
    d <- start(last(rp)) - end(first(rp))
    j <- which(d < 0)
    d[j] <- end(first(rp))[j]-start(last(rp))[j]
    big <- d > 3e3
    thin <- seq(1, length(rp), length.out=max.out)
    thin <- sort(unique(c(thin, which(big))))
    rp <- rp[thin]
  }
  rp
}

meltReadPairs <- function(rps){
  is.improper <- names(rps) != ""
  names(rps) <- NULL
  df <- data.frame(seqnames=as.character(seqnames(rps)),
                  start.first=start(first(rps)), end.first=end(first(rps)),
                  start.last=start(last(rps)),
                  end.last=end(last(rps)))
  ##df <- as(rps, "data.frame")
  df$readpair <- seq_len(nrow(df))
  df1 <- df[, c("start.first", "end.first", "readpair")]
  df2 <- df[, c("start.last", "end.last", "readpair")]
  colnames(df1) <- colnames(df2) <- c("start", "end", "readpair")
  df <- rbind(df1, df2)
  df$read <- rep(c("first", "last"), c(nrow(df1), nrow(df2)))
  df$read <- factor(df$read, levels=c("first", "last"))
  df$is.improper <- rep(is.improper, 2)
  df
}

#' Collect preprocessed bin-level data, segmentation, and improperly paired
#' reads in a list
#'
#' @param bam.file length-one character vector providing path to BAM file
#' @param genome length-one character vector providing genome build (hg18 or hg19)
#' @param bins a \code{GRanges} object with \code{log_ratio} in the \code{mcols}
#' @param segments a \code{GRanges} object with \code{seg.mean} in the \code{mcols}
#' @param improper_rp a \code{GAlignmentPairs} object with improperly paired reads
#' @export
preprocessData <- function(bam.file=NULL,
                           genome,
                           bins,
                           segments,
                           read_pairs){
  list(bam.file=bam.file,
       genome=genome,
       bins=bins,
       segments=segments,
       read_pairs=read_pairs)
}


preprocessData2 <- function(bam.file=NULL,
                           genome,
                           bins,
                           segments,
                           read_pairs){
  list(bam.file=bam.file,
       genome=genome,
       bins=bins,
       segments=segments,
       read_pairs=read_pairs)
}
