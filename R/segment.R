setGeneric("CNAObject", function(object, valuename) standardGeneric("CNAObject"))

setMethod("CNAObject", c("PreprocessViews2", "missing"), function(object, valuename){
  bins <- rowRanges(object)
  x <- assays(object)[, 1, drop=FALSE]
  CNA(x,
      chrom=as.character(seqnames(bins)),
      maploc=start(bins),
      sampleid=colnames(object),
      presorted=TRUE)
})

setMethod("CNAObject", c("GenomicRanges", "character"), function(object, valuename){
  x <- mcols(object)[[valuename]]
  CNA(as.matrix(x),
      chrom=as.character(seqnames(object)),
      maploc=start(object),
      sampleid=object$id[1],
      presorted=TRUE)
})

not_in_filters <- function(x, filters){
  not_filtered <- rep(TRUE, length(x))
  if(!missing(filters)) {
    for(j in seq_along(filters)){
      tmp <- !eval(filters[[j]])
      not_filtered <- tmp & not_filtered
    }
  }
  not_filtered
}


.segment <- function(object, param=SegmentParam(), ...){
    filterList <- list(x=expression(is.na(x)))
    x <- assays(object)[, 1]
    ##files <- file.path(tree[["cbs"]], rdsId(object))
    select <- not_in_filters(x, filterList[["x"]])
    ##x <- x[select]
    dat <- CNAObject(object[select, ])
    seg <- segment(dat,
                   alpha=cbs_alpha(param),
                   undo.splits=undo.splits(param),
                   undo.SD=undo.SD(param),
                   verbose=param@verbose, ...)
    g <- cbs2granges(seg, seqinfo(rowRanges(object)))
    g$sample <- colnames(object)
    g
}

.segmentBins <- function(bins, param, ...){
  ##bins <- bins[!is.na(bins$adjusted)]
  ##adjusted <- bins$adjusted
  bins <- bins[!is.na(bins$log_ratio)]
  adjusted <- bins$log_ratio
  cnaobj <- CNA(as.matrix(adjusted),
                chrom=as.character(seqnames(bins)),
                maploc=start(bins),
                ##sampleid=object$id[1],
                presorted=TRUE)
  seg <- segment(cnaobj,
                 alpha=cbs_alpha(param),
                 undo.splits=undo.splits(param),
                 undo.SD=undo.SD(param),
                 verbose=param@verbose, ...)
  g <- cbs2granges(seg, seqinfo(bins))
  ## g <- tryCatch(cbs2granges(seg, seqinfo(bins)), error=function(e) NULL)
  ## if(is.null(g)) browser()
  g
}

#' Segment log2-transformed and GC-adjusted counts using Circular Binary Segmentation.
#'
#' Segmentation of the transformed counts is performed by the Circular
#' Binary Segmentation (CBS) algorithm implemented in the \code{DNAcopy}
#' package.
#'
#' @param param a \code{SegmentParam} object
#'
#' @param bins a \code{GRanges} object containing adjusted counts. Note, the
#'   adjusted counts must be stored in a column named \code{log_ratio}.
#'
#' @param ... additional arguments to \code{segment} in the \code{DNAcopy} package
#'
#' @details A default set of segmentation parameters are provided in
#' \code{SegmentParam()} and can be altered by creating a new \code{SegmentParam}
#' object if desired.
#'
#' Some users may wish to control additional parameters to the CBS algorithm and
#' therefore \code{segmentBins} is designed to take any argument from \code{DNAcopy::segment}
#' as input (see examples).
#'
#' @seealso \code{\link[DNAcopy]{segment}} for description of circular binary
#'   segmentation and references therein; see
#'   \code{\link[svclasses]{SegmentParam-class}} for a description of the
#'   default parameters settings passed to the \code{segment} function. See
#'   \code{\link[svpreprocess]{binNormalize}} for obtaining normalized counts
#'   for segmentation.
#'
#' @examples
#'   library(svfilters.hg19)
#'   data(bins1kb)
#'   library(GenomeInfoDb)
#'   library(DNAcopy)
#'   bins1kb <- keepSeqlevels(bins1kb, "chr22", pruning.mode = "coarse")
#'   bins1kb$log_ratio <- c(rnorm(ceiling(length(bins1kb)/2), mean = 0, sd = 0.4),
#'                         rnorm(floor(length(bins1kb)/2), mean = -1, sd = 0.4))
#'   segmentBins(bins1kb) # Using default segmentation parameters
#'   segmentBins(bins1kb, param = SegmentParam(alpha = 0.01, undo.splits = "sdundo",
#'                                             undo.SD = 5, verbose = 1))
#'   segmentBins(bins1kb, param = SegmentParam(),
#'               weights = abs(rnorm(length(bins1kb))))
#'
#' @export
segmentBins <- function(bins, param=SegmentParam(), ...){
  gen <- genome(bins)[[1]]
  data(gaps, package=paste("svfilters", gen, sep="."),
       envir=environment())
  gaps <- gaps[gaps$type=="centromere", ]
  seqlevels(gaps, pruning.mode="coarse") <- seqlevels(bins)
  parms <- GRanges(seqnames(gaps), IRanges(rep(1, length(gaps)), start(gaps)))
  parms$arm <- paste0(chromosome(parms), "p")
  qarms <- GRanges(seqnames(gaps), IRanges(end(gaps), seqlengths(gaps)))
  qarms$arm <- paste0(chromosome(qarms), "q")
  arms <- sort(c(parms, qarms))
  bins$arm <- NA
  hits <- findOverlaps(bins, arms)
  bins$arm[queryHits(hits)] <- arms$arm[subjectHits(hits)]
  bins$arm <- factor(bins$arm, levels=arms$arm)
  bins_grl <- split(bins, bins$arm)
  ##  chroms <- seqlevels(bins)
  results <- vector("list", length(bins_grl))
  for(i in seq_along(bins_grl)){
    chrbins <- bins_grl[[i]]
    chrbins <- keepSeqlevels(chrbins, as.character(unique(seqnames(chrbins))))
    if(length(chrbins) < 2) next()
    results[[i]] <- .segmentBins(chrbins, param=param, ...)
  }
  results <- results[ !sapply(results, is.null) ]
  g <- unlist(GRangesList(results))
  seqlevels(g, pruning.mode="coarse") <- seqlevels(bins)
  seqinfo(g) <- seqinfo(bins)
  g
}


segmentCoverage <- function(object, param=SegmentParam(), ...){
  object <- object[, 1]
  chroms <- seqlevels(object)
  results <- vector("list", length(chroms))
  for(i in seq_along(chroms)){
    chr <- chroms[i]
    obj <- keepSeqlevels(object, chr, pruning.mode="coarse")
    if(nrow(obj) < 2) next()
    results[[i]] <- .segment(obj, param=param, ...)
  }
  g <- unlist(GRangesList(results))
  g
}

checkMcols <- function(x, id){
  g <- granges(x)
  g$seg.mean <- x$seg.mean
  g$sample <- rep(id, length(x))
  g
}

#' Segments log2-transformed and GC-adjusted counts.
#'
#' Segmentation of the transformed counts is performed by the circular
#' binary segmentation algorithm implemented in the \code{DNAcopy}
#' package.
#'
#' @seealso \code{\link[DNAcopy]{segment}} for description of circular
#'   binary segmentation and references therein; see
#'   \code{\link[svclasses]{SegmentParam-class}} for a description of the default
#'   parameters settings passed to the \code{segment} function.
#'
#' @return A \code{GRangesList} object.  Each element is the set of
#'   \code{GRanges} for a given sample.  Meta-columns of the
#'   \code{GRanges} elements are \code{seg.mean} (the segment mean)
#'   and \code{sample} (an id for the sample).  Note that
#'   \code{segmentExperiment} saves the \code{GRangesList} object as
#'   an intermediate file in the provided directory tree. Parameters
#'   to the \code{segment} function are passed by an instance of the
#'   \code{SegmentParam} class.
#' 
#' @param object A \code{PreprocessViews2} object
#' @param tree A directory tree for storing intermediate files. See
#'   \code{\link[svclasses]{DataPaths-class}}.
#' @param param An object of class \code{SegmentParam}. 
#' @param ... Additional arguments are passed to the \code{segment}
#'   function in the \code{DNAcopy} package.
segmentExperiment <- function(object, tree, param=SegmentParam(), ...){
  .Deprecated()
  files <- file.path(tree[["cbs"]], rdsId(object))
  J <- seq_len(ncol(object))
  grl <- vector("list", length(J))
  names(grl) <- colnames(object)
  for (j in J){
    file <- files[j]
    if(file.exists(file)){
      g <- readRDS(file)
      grl[[j]] <- g
      next()
    }
    obj <- object[, j]
    g <- segmentCoverage(obj, param=param, ...)
    saveRDS(g, file=file)
    grl[[j]] <- g
  }
  grl2 <- list()
  for(i in seq_along(grl)){
    grl2[[i]] <- checkMcols(grl[[i]], colnames(object)[i])
  }
  grl2 <- GRangesList(grl2)
  names(grl2) <- colnames(object)
  grl2
}
