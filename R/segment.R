constructCNAObject <- function(x, view){
  qr <- rowRanges(view)
  CNA(as.matrix(x),
      chrom=as.character(seqnames(qr)),
      maploc=start(qr),
      sampleid=colnames(view),
      presorted=TRUE)
}

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
    x <- x[select]
    dat <- constructCNAObject(x, object)
    seg <- segment(dat,
                   alpha=cbs_alpha(param),
                   undo.splits=undo.splits(param),
                   undo.SD=undo.SD(param),
                   verbose=param@verbose, ...)
    g <- cbs2granges(seg, seqinfo(rowRanges(object)))
    g$sample <- colnames(object)
    g
}

segmentCoverage <- function(object, param=SegmentParam(), ...){
  object <- object[, 1]
  chroms <- seqlevels(object)
  results <- vector("list", length(chroms))
  for(i in seq_along(chroms)){
    chr <- chroms[i]
    obj <- keepSeqlevels(object, chr)
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
#' @examples
#' library(Rsamtools)
#' library(svpreprocess)
#' extdir <- system.file("extdata", package="svclasses")
#' bviews <- readRDS(file.path(extdir, "bamviews_example.rds"))
#' 
#' set.seed(123)
#' gc <- round(runif(nrow(bviews), 30, 80), 0)
#' mp <- as.integer(round(runif(nrow(bviews), 600, 1000), 0))
#' bamRanges(bviews)$gc <- gc
#' bamRanges(bviews)$map <- mp
#' 
#' unlink(DataPaths(tempdir(), "Test", dryrun=TRUE))
#' tree <- DataPaths(tempdir(), "Test", dryrun=FALSE)
#' views <- countExperiment2(bviews, tree)
#' colnames(views) <- gsub("\\.bam", "", colnames(views))
#' views <- transformExperiment(views, tree)
#' views <- gcExperiment2(views, tree)
#' grl <- segmentExperiment(views, tree)
#' class(grl)
#' unlist(grl)
#' identical(names(grl), colnames(views))
#'
#'
#' ## The granges object for each sample is saved in the provided
#' ##  directory tree.  In particular, list.files(tree[["cbs"]]) ##
#' ##  Calling segmentExperiment a second time merely reads these
#' ##  objects from disk:
#' grl2 <- segmentExperiment(views, tree)
#' identical(grl, grl2)
#'
#' ## To rerun the sementation, these files must be removed or a
#' ##  different 'tree' object must be passed. For example, if we wanted
#' ##  to rerun with different segmentation parameters
#' 
#' unlink(list.files(tree[["cbs"]], full.names=TRUE))
#' sp <- SegmentParam(alpha=0.05)
#' \dontrun{
#'   grl3 <- segmentExperiment(views, tree, sp)
#' }
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
#' @export
#' @param object A \code{PreprocessViews2} object
#' @param tree A directory tree for storing intermediate files. See
#'   \code{\link[svclasses]{DataPaths-class}}.
#' @param param An object of class \code{SegmentParam}. 
#' @param ... Additional arguments are passed to the \code{segment}
#'   function in the \code{DNAcopy} package.
segmentExperiment <- function(object, tree, param=SegmentParam(), ...){
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
