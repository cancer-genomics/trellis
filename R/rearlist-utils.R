#' @include rearrangement-utils.R
NULL

setMethod("numberLinkingRP", "RearrangementList", function(object){
  nrp <- sapply(object, numberLinkingRP)
  as.integer(nrp)
})

#' @aliases linkedTo,GRanges-method
#' @rdname linkedTo
setMethod("linkedTo", "GRanges", function(x)  x$linked.to)

#' @aliases linkedTo,Rearrangement-method
#' @rdname linkedTo
setMethod("linkedTo", "Rearrangement", function(x)  linkedBins(x)$linked.to)

#' @aliases linkedTo,RearrangementList-method
#' @rdname linkedTo
setMethod("linkedTo", "RearrangementList", function(x)  linkedBins(x)$linked.to )


minimumBinSize <- function(object){
  lb <- linkedBins(object)
  bin.sizes <- cbind(width(lb), width(lb$linked.to))
  rowMins(bin.sizes)
}

nearAmplicon <- function(object, amplicons, maxgap=2e3){
  lb <- linkedBins(object)
  h1 <- findOverlaps(amplicons, lb, maxgap=maxgap)
  h2 <- findOverlaps(amplicons, linkedTo(lb), maxgap=maxgap)
  is_amp <- rep(FALSE, length(object))
  if(length(h1) | length(h2) > 0){
    is_amp[c(subjectHits(h1), subjectHits(h2))] <- TRUE
  }
  is_amp
}

nearDeletion <- function(object, del, maxgap=2e3){
  lb <- linkedBins(object)
  h1 <- findOverlaps(del, lb, maxgap=maxgap)
  h2 <- findOverlaps(del, linkedTo(lb), maxgap=maxgap)
  is_del <- rep(FALSE, length(object))
  if(length(h1) | length(h2) > 0){
    is_del[c(subjectHits(h1), subjectHits(h2))] <- TRUE
  }
  is_del
}


#' Applies germline CNV and sequence-based filters for non-redundant
#' identification of germline rearrangements
#'
#' @seealso See \code{\link{filterRearExperiment}} for a wrapper that
#'   saves filtered results. See \code{\link{filterRearrangementList}} for
#'   filtering somatic rearrangements.
#' 
#' @export
#' @param object A \code{RearrangementList} object
#' @param params A \code{RearrangementParams} object
#' @param germline_filters A \code{GRanges} object of germline CNVs and
#'   sequence-based filters
#' @param return_stats length-one logical vector. For internal use only
filterGermlineRear <- function(object, params, germline_filters, return_stats=FALSE){
  stats <- matrix(NA, 4, 2)
  dimnames(stats) <- rev(list(c("n_keep", "n_dropped"),
                              c("fractionLinking",
                                "modalType",
                                "cluster_size",
                                "rearr_filter")))
  if(length(object)==0) {
    stats[,] <- 0
    return(list(rlist=object, stats=stats))
  }
  is_amp <- nearAmplicon(object, GRanges())
  object$isnear_amp <- is_amp
  is_del <- nearDeletion(object, GRanges())
  object$isnear_del <- is_del
  object$prop_modal_call <- percentRearrangement(object)
  drop <-  object$prop_modal_call < percentModalType(params)
  idrop <- as.integer(table(factor(drop, levels=c("FALSE", "TRUE"))))
  stats[2, ] <- idrop
  if(any(drop))  object <- object[ !drop ]
  if(length(object) == 0){
    if(return_stats) return(stats)
    return(object)
  }
  object$min_bin_size <- minimumBinSize(object)
  drop <-   object$min_bin_size < minClusterSize(params)
  idrop <- as.integer(table(factor(drop, levels=c("FALSE", "TRUE"))))
  stats[3, ] <- idrop
  if(any(drop)) object <- object[ !drop ]
  if(length(germline_filters) > 0){
    lb <- linkedBins(object)
    overlaps_other1 <- overlapsAny(lb, germline_filters)
    overlaps_other2 <- overlapsAny(lb$linked.to, germline_filters)
    overlaps_other <- overlaps_other1 | overlaps_other2
    object$overlaps_cnv_filter <- overlaps_other
  }
  object$number_linking <- numberLinkingRP(object)
  if(return_stats) return(stats)
  object

}

#' Applies several filters for somatic rearrangements
#' 
#' Removes rearrangement intervals that overlap with previously
#' identified somatic amplicons and deletions, as well as germline
#' CNVs and rearrangements.
#'
#'  
#' @export
#' @param object A \code{RearrangementList}
#' @param params A \code{RearrangementParams} object
#' @param amps \code{GRanges} of amplicons
#' @param del.gr \code{GRanges} of deletions
#' @param rear_filter \code{GRanges} object of germline rearrangements
#' @param germline_filters \code{GRanges} of germline CNVs and sequence filters
#' @param return_stats a length-one logical vector.  For internal use.
filterRearrangementList <- function(object, params,
                                    amps, del.gr,
                                    rear_filter,
                                    germline_filters=GRanges(),
                                    return_stats=FALSE){
  stats <- matrix(NA, 4, 2)
  dimnames(stats) <- rev(list(c("n_keep", "n_dropped"),
                              c("fractionLinking",
                                "modalType",
                                "cluster_size",
                                "rearr_filter")))

  if(length(object)==0) {
    stats[,] <- 0
    x <- list(rlist=object, stats=stats)
    return(x)
  }
  ## for each element, check for overlap with previously identified
  ## deletions and amplifications
  if(!missing(amps)){
    is_amp <- nearAmplicon(object, amps)
    object$isnear_amp <- is_amp
  }
  if(!missing(del.gr)){
    is_del <- nearDeletion(object, del.gr)
    object$isnear_del <- is_del
    ##
    ## Shouldn't this be !is_amp instead of is_del
    ## - is_del is just an indicator for whether the junctionn is near
    ## the deletion GRanges filter (del.gr)
    drop <- is_del & fractionLinkingTags(object) < percentLinking(params)
  } else{
    ## this might drop amplifications that are real
    drop <- fractionLinkingTags(object) < percentLinking(params)
  }
  idrop <- as.integer(table(factor(drop, levels=c("FALSE", "TRUE"))))
  stats[1, ] <- idrop
  if(any(drop))  object <- object[ !drop  ]

  object$prop_modal_call <- percentRearrangement(object)
  drop <-  object$prop_modal_call < percentModalType(params)
  idrop <- as.integer(table(factor(drop, levels=c("FALSE", "TRUE"))))
  stats[2, ] <- idrop
  if(any(drop))  object <- object[ !drop ]

  if(length(object) == 0){
    if(return_stats) return(stats)
    return(object)
  }
  object$min_bin_size <- minimumBinSize(object)
  drop <-   object$min_bin_size < minClusterSize(params)
  idrop <- as.integer(table(factor(drop, levels=c("FALSE", "TRUE"))))
  stats[3, ] <- idrop
  if(any(drop)) object <- object[ !drop ]
  ## flag any overlap with rearrangement filter

  ##object$overlaps_rearrangement
  if(!missing(rear_filter)){
    drop <- overlapsAny(object, rear_filter)
    idrop <- as.integer(table(factor(drop, levels=c("FALSE", "TRUE"))))
    stats[4, ] <- idrop
    if(any(drop)) object <- object[ !drop ]
  }
  if(length(object)==0){
    if(return_stats) return(stats)
    return(object)
  }
  if(length(germline_filters) > 0){
    lb <- linkedBins(object)
    overlaps_other1 <- overlapsAny(lb, germline_filters)
    overlaps_other2 <- overlapsAny(lb$linked.to, germline_filters)
    overlaps_other <- overlaps_other1 | overlaps_other2
    object$overlaps_cnv_filter <- overlaps_other
  }
  object$number_linking <- numberLinkingRP(object)
  ##list(rlist=object, stats=stats)
  if(return_stats) return(stats)
  object
}




filterRearrangementList2 <- function(rlist,
                                     params=RearrangementParams(),
                                     amps,
                                     del.gr,
                                     rear_filter,
                                     germline_filters=GRanges(),
                                     return_stats=FALSE){
  stats <- matrix(NA, 4, 2)
  dimnames(stats) <- rev(list(c("n_keep", "n_dropped"),
                              c("fractionLinking",
                                "modalType",
                                "cluster_size",
                                "rearr_filter")))

  if(length(rlist)==0) {
    stats[,] <- 0
    x <- list(rlist=rlist, stats=stats)
    return(x)
  }
  ## for each element, check for overlap with previously identified
  ## deletions and amplifications
  if(!missing(amps)){
    rlist$isnear_amp <- nearAmplicon(rlist, amps)
  }
  if(!missing(del.gr)){
    is_del <- nearDeletion(rlist, del.gr)
    rlist$isnear_del <- is_del
    ##
    ## Shouldn't this be !is_amp instead of is_del
    ## - is_del is just an indicator for whether the junctionn is near
    ## the deletion GRanges filter (del.gr)
    drop <- is_del & fractionLinkingTags(rlist) < percentLinking(params)
  } else{
    ## this might drop amplifications that are real
    drop <- fractionLinkingTags(rlist) < percentLinking(params)
  }
  idrop <- as.integer(table(factor(drop, levels=c("FALSE", "TRUE"))))
  stats[1, ] <- idrop
  if(any(drop))  rlist <- rlist[ !drop  ]

  rlist$prop_modal_call <- percentRearrangement(rlist)
  drop <-  rlist$prop_modal_call < percentModalType(params)
  idrop <- as.integer(table(factor(drop, levels=c("FALSE", "TRUE"))))
  stats[2, ] <- idrop
  if(any(drop))  rlist <- rlist[ !drop ]

  if(length(rlist) == 0){
    if(return_stats) return(stats)
    return(rlist)
  }
  rlist$min_bin_size <- minimumBinSize(rlist)
  drop <-   rlist$min_bin_size < minClusterSize(params)
  idrop <- as.integer(table(factor(drop, levels=c("FALSE", "TRUE"))))
  stats[3, ] <- idrop
  if(any(drop)) rlist <- rlist[ !drop ]
  ## flag any overlap with rearrangement filter

  ##rlist$overlaps_rearrangement
  if(!missing(rear_filter)){
    drop <- overlapsAny(rlist, rear_filter)
    idrop <- as.integer(table(factor(drop, levels=c("FALSE", "TRUE"))))
    stats[4, ] <- idrop
    if(any(drop)) rlist <- rlist[ !drop ]
  }
  if(length(rlist)==0){
    if(return_stats) return(stats)
    return(rlist)
  }
  if(length(germline_filters) > 0){
    lb <- linkedBins(rlist)
    overlaps_other1 <- overlapsAny(lb, germline_filters)
    overlaps_other2 <- overlapsAny(lb$linked.to, germline_filters)
    overlaps_other <- overlaps_other1 | overlaps_other2
    rlist$overlaps_cnv_filter <- overlaps_other
  }
  rlist$number_linking <- numberLinkingRP(rlist)
  ##list(rlist=rlist, stats=stats)
  if(return_stats) return(stats)
  rlist
}

dropRearNearDeletion <- function(rlist, del.gr, params){
  if(length(del.gr) > 0){
    is_del <- nearDeletion(rlist, del.gr)
    rlist$isnear_del <- is_del
    ##
    ## Shouldn't this be !is_amp instead of is_del
    ## - is_del is just an indicator for whether the junctionn is near
    ## the deletion GRanges filter (del.gr)
    drop <- is_del & fractionLinkingTags(rlist) < percentLinking(params)
  } else{
    ## this might drop amplifications that are real
    drop <- fractionLinkingTags(rlist) < percentLinking(params)
  }
  if(any(drop))  rlist <- rlist[ !drop  ]
  rlist
}

overlapsCNV <- function(rlist, filters, params){
  if(length(filters$amplicons) > 0) rlist$isnear_amp <- nearAmplicon(rlist, filters$amplicons)
  rlist <- dropRearNearDeletion(rlist, filters$deletions, params)
  rlist
}

rearStats <- function(rlist){
  rlist$prop_modal_call <- percentRearrangement(rlist)
  rlist$min_bin_size <- minimumBinSize(rlist)
  rlist
}

applyFilters <- function(rlist, filters, params){
  rlist <- overlapsCNV(rlist, filters, params)
  rlist <- rearStats(rlist)
  is.consistent <- rlist$prop_modal_call >= percentModalType(params)
  is.minsize <- rlist$min_bin_size >= minClusterSize(params)
  overlaps.rear <- overlapsAny(rlist, filters$rear)
  rlist$overlaps_cnv_filter <- overlapsGermlineRear(rlist, filters$germline)
  rlist[is.consistent & is.minsize & !overlaps.rear]
}

#' Filter candidate rearrangements for germline
#'
#' @param rdata a list of rearrangements and filters
#' @param params a \code{RearrangementParams} object
#' @export
filterRear <- function(rdata, params=RearrangementParams()){
  rlist <- rdata$rlist
  if(length(rlist)==0) return(rlist)
  rlist <- applyFilters(rdata$rlist, rdata$filters, params)
  if(length(rlist) == 0) return(rlist)
  ## flag any overlap with rearrangement filter
  rlist$number_linking <- numberLinkingRP(rlist)
  ## for compatibility with filterRearrangementList2
  rlist@modal_rearrangement <- character()
  rlist@percent_rearrangement <- numeric()
  rlist
}

overlapsGermlineRear <- function(rlist, germline_filters){
  lb <- linkedBins(rlist)
  overlaps_other1 <- overlapsAny(lb, germline_filters)
  overlaps_other2 <- overlapsAny(lb$linked.to, germline_filters)
  overlaps_other1 | overlaps_other2
}


#' Exclude linked tag cluster that overlap with germline CNVs or other
#' sequence artifacts
#'
#' This is a wrapper for \code{filterRearrangementList} (somatic) and
#' \code{filterGermlineRear} (germline).  The filtered rearrangements
#' are saved to an intermediate file.  If the file already exists, the
#' saved result is read from disk.  To rerun the filtering step, the
#' intermediate files must be removed.  See examples.
#'
#' REFACTORING: Too many arguments. Unclear distinction between
#' germline_filters and rear_filter.  Should have a separate function for
#' germline, or have an S4 germline class.
#'
#' @seealso See \code{filterRearrangementList} for filtering somatic
#'   rearrangements and \code{filterGermlineRear} for filtering
#'   germline rearrangements.  See \code{\link{listCNVs}} for listing
#'   the deletions and amplicons for a given sample.
#'
#' @details
#' 
#' The amplicons (amps) and deletions (del.gr) will be subset by the
#'   mcols variable 'id'.
#'
#' @return a \code{RearrangementList}
#' @export
#' @param x a \code{RearrangementList}
#' @param dirs a \code{DataPaths} object
#' @param cnvs a list of deletions and amplicons
#' @param germline_filters a \code{GRanges} object of germline filters
#' @param rear_filter a \code{GRanges} object of rearrangements found in germline
#' @param rp a \code{RearrangementParams} object
filterRearExperiment <- function(x,  
                                 dirs,
                                 cnvs,
                                 germline_filters,
                                 rear_filter,
                                 rp=RearrangementParams()){
  amps <- cnvs[["amplicons"]]
  dels <- cnvs[["deletions"]]
  id <- names(x)
  ffiles <- file.path(dirs[["1somatic"]], paste0(id, ".rds"))
  rlist <- setNames(vector("list", length(x)), id)
  names(ffiles) <- id
  for(i in seq_along(id)){
    idx <- id[i]
    outfile <- ffiles[idx]
    if(file.exists(outfile)) {
      rlist[[i]] <- readRDS(outfile)
      next()
    }
    tmp <- filterRearrangementList(x[[i]],
                                   params=rp,
                                   amps=amps,
                                   del.gr=dels,
                                   rear_filter=rear_filter,
                                   germline_filters=germline_filters)
    rlist[[i]] <- tmp
    saveRDS(tmp, file=outfile)
  }
  rlist
}

#' Exclude linked tag cluster that overlap with germline CNVs or other
#' sequence artifacts
#' 
#'  This function is intended for internal use, primarily for
#'  delineating rearrangements that may occur in the germline.
#'
#' @keywords internal
#' @export
#' @param x a \code{RearrangementList}
#' @param dirs a \code{DataPaths} object
#' @param rp a \code{RearrangementParams} object
#' @param germline_filters a \code{GRanges} object of germline CNV and sequence filters
filterGermlineRearList <- function(x,
                                   dirs,
                                   rp=RearrangementParams(),
                                   germline_filters){
  id <- names(x)
  ffiles <- file.path(dirs[["rearrangements/germline"]], paste0(id, ".rds"))
  rlist <- setNames(vector("list", length(x)), id)
  names(ffiles) <- id
  for(i in seq_along(id)){
    idx <- id[i]
    outfile <- ffiles[idx]
    if(file.exists(outfile)) {
      rlist[[i]] <- readRDS(outfile)
      next()
    }
    tmp <- filterGermlineRear(x[[i]],
                              params=rp,
                              germline_filters=germline_filters)
    rlist[[i]] <- tmp
    saveRDS(tmp, file=outfile)
  }
  rlist
}

#' Construct list of filters for somatic rearrangement analyses
#'
#' @param germline a \code{GRanges} object of germline filters
#' @param rear a \code{GRanges} object of rearrangements identified in germline samples (e.g., possibly artifactual)
#' @param deletions a \code{GRanges} object of deletions
#' @param amplicons a \code{GRanges} object of amplicons
#' @return a \code{GRangesList} of rearrangement filter
#' @export
rFilters <- function(germline=GRanges(),
                     rear=GRanges(),
                     deletions=GRanges(),
                     amplicons=GRanges()){
  GRangesList(germline=granges(germline),
              rear=granges(rear),
              deletions=granges(deletions),
              amplicons=granges(amplicons))
}

#' Creates a list object of data required to identify rearrangements
#'
#' @param bins a \code{GRanges} object of non-overlapping (typically 1kb) bins
#' @param read_pairs a \code{GAlignmentPairs}-representation of paired reads
#' @param rlist a \code{RearrangementList}
#' @param filters a \code{GRangesList} of sequence- and germline-based filters
#' @return a list
#' @export
rearrangementData <- function(bins=GRanges(),
                              read_pairs,
                              rlist=RearrangementList(),
                              filters=rFilters()){
  list(bins=bins,
       rlist=rlist,
       read_pairs=read_pairs,
       filters=filters)
}


annotateSide <- function(alignments, regions){
  regions <- sort(regions)
  i <- findOverlaps(alignments, regions, select="first")
  c("left", "right")[i]
}

clusterRegions <- function(rear){
  c(granges(linkedBins(rear)), linkedTo(rear))
}

improperSides <- function(rear){
  irp <- improper(rear)
  f <- first(irp)
  l <- last(irp)
  regions <- clusterRegions(rear)
  mcols(f)$side <- annotateSide(f, regions)
  mcols(l)$side <- annotateSide(l, regions)
  first(irp) <- f
  last(irp) <- l
  ##rear@improper <- irp
  ##rear
  irp
}

splitreadSide <- function(split_reads, regions){
  split_reads$side <- annotateSide(split_reads, regions)
  split_reads
}

.organize_reads <- function(rear, split_reads){
  galp <- improperSides(rear)
  clust1 <- as(first(galp), "GRanges")
  clust1$read <- "R1"
  clust2 <- as(last(galp), "GRanges")
  clust2$read <- "R2"
  names(clust1) <- names(clust2) <- NULL
  ## read pair id
  rpid <- rep(seq_along(clust1), 2)
  clust.gr <- c(clust1, clust2)
  sides <- clust.gr$side
  read <- clust.gr$read
  clust.gr <- granges(clust.gr)
  clust.gr$side <- sides
  clust.gr$rpid <- rpid
  clust.gr$read <- read
  clust.gr$type <- "improper"
  if(!missing(split_reads)){
    sr.side <- annotateSide(split_reads, clusterRegions(rear))
    sr.rpid <- as.numeric(factor(split_reads$qname)) + max(rpid)
    sr.gr <- granges(split_reads)
    sr.gr$side <- sr.side
    sr.gr$rpid <- sr.rpid
    sr.gr$read <- "split"
    sr.gr$type <- "split read"
    gr <- c(clust.gr, sr.gr)
  } else gr <- clust.gr
  df <- as(gr, "data.frame")
  df$side <- factor(df$side)
  ##df$rpid <- factor(df$rpid)
  df
}

.organize_with_split_reads <- function(rlist, split_reads){
  results <- vector("list", length(rlist))
  split_reads <- split_reads[names(rlist)]
  for(i in seq_along(rlist)){
    tmp <- .organize_reads(rlist[[i]], split_reads[[i]])
    tmp$rid <- names(rlist)[i]
    results[[i]] <- tmp
  }
  results
}

.organize_without_split_reads <- function(rlist){
  results <- vector("list", length(rlist))
  for(i in seq_along(rlist)){
    tmp <- .organize_reads(rlist[[i]])
    tmp$rid <- names(rlist)[i]
    results[[i]] <- tmp
  }
  results
}

#' Organize improper read pairs supporting a rearrangement into a data.frame
#'
#' Facilitates plotting reads supporting a rearrangement with ggplot.
#'
#' @param rlist a `RearrangementList`
#' @param split_reads a data.frame of split reads from BLAT. See `rearrangedReads`.
#' @seealso \code{\link{rearrangedReads}}
#' @export
organizeReads <- function(rlist, split_reads){
  if(!missing(split_reads)){
    results <- .organize_with_split_reads(rlist, split_reads)
    return(results)
  }
  .organize_without_split_reads(rlist)
}

bp_breaks <- function(n=5, unit=1000){
  function(x) {
    rng <- 1/unit * range(x, na.rm = TRUE)
    min <- floor(rng[1])
    max <- ceiling(rng[2])
    if (max == min)
      return(unit*min)
    by <- floor((max - min)/n) + 1
    unit*seq(min, max, by = by)
  }
}

bp_trans <- function(unit){
  if(unit == 1000){
    symb <- "kb"
  }
  if(unit == 1e6){
    symb <- "Mb"
  }
  trans <- function(x) x/unit
  inv <- function(x) x*unit
  trans_new(symb, trans, inv,
            bp_breaks(unit = unit),
            domain = c(1e-100, Inf))
}

kb_trans <- function() bp_trans(1000)

mb_trans <- function() bp_trans(1e6)

scale_x_kb <- function(...){
  scale_x_continuous(..., trans=kb_trans())
}

#' Use MB scale in gg-style plots
#'
#' @param ... additional arguments to scale_x_continuous
#' @export
scale_x_mb <- function(...){
  scale_x_continuous(..., trans=mb_trans())
}

readOrientation <- function(){
  ## same chromosome
  ## R1+ < R2-:  stays the same
  ## R1+ > R2-:  flip orientation
  ## R1- > R2+:  stays the same
  ##
  ## different chromosome
  ## R1+,R2-:  
  ##
}

fiveTo3PrimeList <- function(rlist, build, maxgap=5000){
  rlist2 <- lapply(rlist, fiveTo3Prime, build=build, maxgap=5000)
  rlist2 <- rlist2[ !sapply(rlist2, is.null) ]
  rlist2
}

#' Loads TxDb object from TxDb.Hsapiens.UCSC.<build>.refGene
#'
#' @return A \code{TxDb} object
#' @param build genome build
#' @export
#' @examples
#' txdb <- loadTxDb("hg19")
#' txdb
loadTxDb <- function(build){
  txdb.pkg <-   paste0("TxDb.Hsapiens.UCSC.", build, ".refGene")
  txdb.path <- system.file("extdata", package=txdb.pkg)
  txdb.file <- file.path(txdb.path, paste0(txdb.pkg, ".sqlite"))
  txdb <- loadDb(txdb.file, packageName=txdb.pkg)
  conn <- dbFileConnect(txdb.file)
  dbFileDisconnect(conn)
  txdb
}

#' Load transcripts, CDS, and a txdb object for a user-specified build
#'
#' @export
#' @param build character string providing genome build. Currently only hg19 and hg18 are supported.
#' @param seq_levels allowed seqnames
#' @return a named list with elements `transcripts` (\code{GRanges} of transcripts), `cds` (\code{GRangesList} of CDS listed by transcript), and `txdb` (a \code{txdb} object).
#' @examples
#' tx.cds <- loadTxdbTranscripts("hg19")
#' tx.cds[["cds"]]
#' tx.cds[["transcripts"]]
#' tx.cds[["txdb"]]
loadTxdbTranscripts <- function(build, seq_levels){
  txdb <- loadTxDb(build)
  if(!missing(seq_levels)){
    txdb <- keepSeqlevels(txdb, seq_levels)
  }
  tx <- transcripts(txdb)
  tx <- tx[grep("NM_", tx$tx_name)]
  cds <- suppressWarnings(cdsBy(txdb, "tx", use.names=TRUE))
  list(transcripts=tx, cds=cds, txdb=txdb)
}

#' Evaluates whether sequence junction of a rearrangement is within 5000bp of a transcript
#'
#' Only rearrangements where both sides of the sequence junction are near a transcript are evaluated as candidate in-frame fusions.
#'
#' @param rlist a \code{RearrangementList}
#' @param build genome build
#' @param maxgap integer for number of basepairs allowed between a candidate sequence junction and the improperly paired read alignments
#' @export
#' @return a logical vector of the same length as the \code{rlist} object
#' @examples
#'  extdata <- system.file("extdata", package="svbams")
#'  rfile <- file.path(extdata, "CGOV11T_1.bam.rds")
#'  rlist <- readRDS(rfile)
#'  head(seqJunctionNearTx(rlist, "hg19"))
seqJunctionNearTx <- function(rlist, build, maxgap=5000){
  tx.cds <- loadTxdbTranscripts(build, seq_levels=seqlevels(rlist))
  tx <- tx.cds[["transcripts"]]
  lb <- linkedBins(rlist)
  region1 <- lb
  region2 <- linkedTo(lb)
  cds.region1 <- overlapsAny(region1, tx, ignore.strand=TRUE, maxgap=maxgap)
  cds.region2 <- overlapsAny(region2, tx, ignore.strand=TRUE, maxgap=maxgap)
  cds.region1 & cds.region2
}

#' Reorders each rearrangement as two 5-prime to 3-prime orientations
#'
#' @param rlist a \code{RearrangementList}
#' @param build character string providing genome build (only hg19 or hg18 are supported)
#' @return a \code{RearrangementList}
#' @export
#' @examples
#'  library(trellis)
#'  extdata <- system.file("extdata", package="svbams")
#'  rfile <- file.path(extdata, "CGOV11T_1.bam.rds")
#'  rlist <- readRDS(rfile)
#'  near.coding <- seqJunctionNearTx(rlist, "hg19")
#'  \dontrun{
#'    rlist2 <- fiveTo3List(rlist[near.coding], build="hg19")
#'    rlist2
#' }
fiveTo3List <- function(rlist, build){
  rlist2 <- RearrangementList()
  nms.list <- vector("list", length(rlist) * 2)
  k <- 1
  for(i in seq_along(rlist)){
    r5p3 <- fiveTo3Prime(rlist[[i]], build)
    rlist2[[k]] <- r5p3[[1]]
    rlist2[[k+1]] <- r5p3[[2]]
    nms.list[[i]] <- names(r5p3)
    k <- k+2
  }
  nms <- unlist(nms.list)
  rlist2@names <- nms
  ##
  ## duplicate colData
  ##
  cd <- rlist@colData
  index <- rep(seq_len(nrow(cd)), each=2)
  cd <- cd[index, ]
  rownames(cd) <- nms
  rlist2@colData <- cd
  ##
  ## remove rearrangements that are not near a gene (why was this not already caught?)
  ##
  lb <- tryCatch(linkedBins(rlist2), error=function(e) NULL)
  if(is.null(lb)){
    lb.list <- lapply(rlist2, linkedBins)
    nc <- sapply(lb.list, function(x) ncol(mcols(x)))
    j <- which(nc != 3)
    ## we didn't determine a gene_name or whether orientation is reversed for these rearrangements
    lb.list2 <- lapply(lb.list[j], function(x){
      mc <- mcols(x)
      mc$gene_name <- as.character(NA)
      mc$reverse <- as.logical(NA)
      mc <- mc[c("gene_name", "reverse", "linked.to")]
      mcols(x) <- mc
      x
    })
    rlist3 <- rlist2[j]
    for(i in seq_along(rlist3)){
      linkedBins(rlist3@data[[i]]) <- lb.list[[i]]
      k <- j[i]
      rlist2[[k]] <- rlist3[[i]]
    }
  }
  rlist2 <- rlist2[ !is.na(linkedBins(rlist2)$gene_name) ]
  rlist3 <- rlist2[ elementNROWS(splitReads(rlist2)) > 0 ]
  nna <- sapply(splitReads(rlist3), function(x) all(is.na(mcols(x)$reverse)))
  rlist4 <- rlist3[ !nna ]
  rlist4
}

numberSplitReads <- function(rlist){
  N <- rep(NA, length(rlist))
  names(N) <- names(rlist)
  for(i in seq_along(rlist)){
    r <- rlist[[i]]
    bins <- linkedBins(r)
    sr1 <- countOverlaps(splitReads(r), bins, maxgap=50)
    sr2 <- countOverlaps(splitReads(r), bins$linked.to, maxgap=50)
    N[i] <- min(sum(sr1), sum(sr2))
  }
  N
}

seqJunction <- function(r, maxgap=50){
  bins <- linkedBins(r)
  sr1 <- subsetByOverlaps(splitReads(r), bins, maxgap=maxgap)
  if(sr1$reverse[1]){
    ## minus strand -- end is actually the start
    ends <- start(sr1)
    s1 <- "-"
  } else{
    ends <- end(sr1)
    s1 <- "+"
  }
  jxn1 <- median(ends)
  sr2 <- subsetByOverlaps(splitReads(r), bins$linked.to, maxgap=maxgap)
  if(sr2$reverse[1]){
    ## minus strand -- start is actually the end
    starts <- end(sr2)
    s2 <- "-"
  } else {
    starts <- start(sr2)
    s2 <- "+"
  }
  jxn2 <- median(starts)
  g <- c(GRanges(chromosome(sr1)[1], IRanges(jxn1, width=1),
                 strand=s1,
                 seqinfo=seqinfo(bins)),
         GRanges(chromosome(sr2)[1], IRanges(jxn2, width=1),
                 strand=s2,
                 seqinfo=seqinfo(bins)))
  names(g) <- c("5p", "3p")
  g
}

#' Use split reads for each rearrangement to more precisely define sequence boundry
#'
#' This function constructs a \code{GRanges} object of the 5-prime and 3-prime sequence junctions for a \code{RearrangementList} that has already been ordered by its two 5-prime to 3-prime orientation.
#'
#' @param rlist a \code{RearrangementList}
#' @return a \code{GRanges} object of the 5-prime and 3-prime sequence junctions.  The 3-prime junctions are included in the '3p' field of the \code{GRanges} object
#' @export
#' @examples
#'  library(trellis)
#'  extdata <- system.file("extdata", package="svbams")
#'  rfile <- file.path(extdata, "CGOV11T_1.bam.rds")
#'  rlist <- readRDS(rfile)
#'  near.coding <- seqJunctionNearTx(rlist, "hg19")
#'  ## Just do the first two
#'  index <- which(near.coding)[1:2]
#'  rlist2 <- fiveTo3List(rlist[index], build="hg19")
#'  jxns <- seqJunctions_Rlist(rlist2)
#'  jxns
seqJunctions_Rlist <- function(rlist){
  jxn.grl <- GRangesList(lapply(rlist, seqJunction))
  names(jxn.grl) <- names(rlist)
  jxn.gr <- unlist(jxn.grl, use.names=FALSE)
  jxn.gr$rid <- rep(names(jxn.grl), each=2)
  index.5p <- seq(1, length(jxn.gr), 2)
  jxn.5p <- jxn.gr[index.5p]
  jxn.3p <- jxn.gr[-index.5p]
  jxn.5p$"3p" <- jxn.3p
  jxn.5p
}

.invalid_split_reads <- function(r, maxgap=100){
  is_invalid <- FALSE
  bins <- linkedBins(r)
  sr1 <- subsetByOverlaps(splitReads(r), bins, maxgap=maxgap)
  sr2 <- subsetByOverlaps(splitReads(r), bins$linked.to, maxgap=maxgap)
  ## no split reads
  if(length(sr1) == 0 && length(sr2) == 0) return(FALSE)
  if(length(sr1) == 0 || length(sr2) == 0){
    is_invalid <- TRUE
  }
  is_invalid
}



is_invalid_splits <- function(rlist, maxgap=100){
  sapply(rlist, .invalid_split_reads, maxgap=maxgap)
}

#' Check that the split reads available for each rearrangement span both sides of the sequence junction
#'
#' @param rlist a \code{RearrangementList}
#' @param maxgap length-one numeric vector passed to \code{subsetByOverlaps}
#' @return logical vector of the same length as \code{rlist}
#' @export
is_valid_splits <- function(rlist, maxgap=100) !is_invalid_splits(rlist, maxgap)
