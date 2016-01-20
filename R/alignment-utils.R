## Adapted from makeGAlignmentPairs in GenomicAlignments
makeGAlignmentPairs2 <- function(x, use.names=FALSE, use.mcols=FALSE, strandMode=1){
  if (!isTRUEorFALSE(use.names))
    stop("'use.names' must be TRUE or FALSE")
  if (!isTRUEorFALSE(use.mcols)) {
    if (!is.character(use.mcols))
      stop("'use.mcols' must be TRUE or FALSE or a character vector ",
           "specifying the metadata columns to propagate")
    if (!all(use.mcols %in% colnames(mcols(x))))
      stop("'use.mcols' must be a subset of 'colnames(mcols(x))'")
  }
  mate <- findMateAlignment(x)
  x_is_first <- GenomicAlignments:::.isFirstSegment.GAlignments(x)
  x_is_last <- GenomicAlignments:::.isLastSegment.GAlignments(x)
  first_idx <- which(!is.na(mate) & x_is_first)
  last_idx <- mate[first_idx]
  ## Fundamental property of the 'mate' vector: it's a permutation of order
  ## 2 and with no fixed point on the set of indices for which 'mate' is
  ## not NA.
  ## Check there are no fixed points.
  if (!all(first_idx != last_idx))
    stop("findMateAlignment() returned an invalid 'mate' vector")
  ## Check order 2 (i.e. permuting a 2nd time brings back the original
  ## set of indices).
  if (!identical(mate[last_idx], first_idx))
    stop("findMateAlignment() returned an invalid 'mate' vector")
  ## One more sanity check.
  if (!all(x_is_last[last_idx]))
    stop("findMateAlignment() returned an invalid 'mate' vector")
  ## Check the 0x2 bit (isProperPair).
  x_is_proper <- as.logical(bamFlagAsBitMatrix(mcols(x)$flag,
                                               bitnames="isProperPair"))
  ans_is_proper <- x_is_proper[first_idx]
  ans_first <- x[first_idx]
  ans_last <- x[last_idx]
  ans_names <- NULL
  if (use.names)
    ans_names <- names(ans_first)
  names(ans_first) <- names(ans_last) <- NULL
  if (is.character(use.mcols)) {
    mcols(ans_first) <- mcols(ans_first)[use.mcols]
    mcols(ans_last) <- mcols(ans_last)[use.mcols]
  } else if (!use.mcols) {
    mcols(ans_first) <- mcols(ans_last) <- NULL
  }
  gps <- GAlignmentPairs(first=ans_first,
                         last=ans_last,
                         isProperPair=ans_is_proper,
                         names=ans_names,
                         strandMode=strandMode)
  gps
}

.trimInvalidReadsGAlign <- function(x){
  chromosome <- function(x) as.character(seqnames(x))
  ends <- setNames(end(x), chromosome(x))
  ends.info <- seqlengths(x)[names(ends)]
  is_valid <- ends <= ends.info
  x[is_valid]
}

#' Creates a default set of flags for reading improperly paired alignments
#'
#' This function is a wrapper to \code{scanBamFlag} in \code{Rsamtools}.
#'
#' @seealso \code{\link[Rsamtools]{scanBamFlag}} for complete details
#'   and \code{\link{improperAlignmentParams}} for a wrapper of this
#'   function that generates a \code{ScanBamParam} object using these
#'   flags.
#' 
#' @examples
#' require(Rsamtools)
#' flags <- scanBamFlag(isDuplicate=FALSE,
#'                      isPaired=TRUE,
#'                      isUnmappedQuery=FALSE,
#'                      hasUnmappedMate=FALSE,
#'                      isProperPair=FALSE)
#' flags2 <- improperAlignmentFlags()
#' identical(flags, flags2)
#' print(flags)
#' 
#' @export
improperAlignmentFlags <- function(){
  flags <- scanBamFlag(isDuplicate=FALSE,
                       isPaired=TRUE,
                       isUnmappedQuery=FALSE,
                       hasUnmappedMate=FALSE,
                       isProperPair=FALSE)
}


#' Helper function for specifying flags for reading improperly
#' paired reads
#'
#' This is a wrapper to \code{ScanBamParam}
#'
#' @seealso \code{\link[Rsamtools]{ScanBamParam}} and
#'   \code{\link[Rsamtools]{bamFlag}} for full documentation in
#'   \code{Rsamtools}.  See \code{improperAlignmentFlags} for the
#'   default set of flags.
#'
#' @examples
#' require(Rsamtools)
#' flags <- improperAlignmentFlags()
#' print(flags)
#' params <- ScanBamParam(flag = flags, what=c("flag", "mrnm", "mpos", "mapq"))
#' params2 <- improperAlignmentParams()
#' print(params2)
#' identical(params, params2)
#' 
#' @param flags A length-two integer vector as provided by \code{improperAlignmentFlags}
#' @export
improperAlignmentParams <- function(flags=improperAlignmentFlags()){
  ScanBamParam(flag=improperAlignmentFlags(), what=c("flag", "mrnm", "mpos", "mapq"))
}

getImproperAlignmentPairs <- function(object,
                                      param,
                                      mapq_thr=-Inf,
                                      use.mcols=TRUE){
  bam.file <- bamPaths(object)
  out.file <- improperPaths(object)
  flags <- improperAlignmentFlags()
  if(missing(param)){
    param <- ScanBamParam(flag=flags, what=c("flag", "mrnm", "mpos", "mapq"))
  }
  irp <- readGAlignments(bam.file, use.names=TRUE, param=param)
  irp <- .trimInvalidReadsGAlign(irp)
  if(mapq_thr > -Inf){
    irp <- irp[ mcols(irp)$mapq >= mapq_thr ]
  }
  irp2 <- makeGAlignmentPairs2(irp, use.mcols=use.mcols, use.names=TRUE)
  saveRDS(irp2, file=out.file)
  irp2
}

#' Parse improper read pairs from a bam file
#'
#' Parses improper read pairs from a bam file and saves the result as
#' a serialized R object.  The file paths to the improper read pairs
#' are given by accessors defined from the \code{AlignmentViews2}
#' class.
#'
#' @examples
#'   library(Rsamtools)
#'   require(TestBams)
#'   extdata <- system.file("extdata", package="TestBams")
#'   bam.file <- list.files(extdata, pattern="\\.bam$", full.name=TRUE)
#'   bv <- BamViews(bam.file)
#'   dp <- DataPaths(tempdir())
#'   aviews <- AlignmentViews2(bv, dp)
#'   \dontrun{
#'     writeImproperAlignments2(aviews)
#'   }
#' 
#' @rdname AlignmentViews2
#' @export
#' @param aviews a \code{AlignmentViews2} object
#' @param param a \code{ScanBamParam} object
#' @param mapq_thr  length-one numeric vector indicating lower limit of MAPQ score
#' @param use.mcols length-one logical vector
writeImproperAlignments2 <- function(aviews,
                                     param=improperAlignmentParams(),
                                     mapq_thr=-Inf, use.mcols=TRUE){
  if(all(file.exists(improperPaths(aviews)))){
    return(invisible())
  }
  aviews <- aviews[, !file.exists(improperPaths(aviews))]
  if(ncol(aviews)==0) return(invisible())
  aln_paths <- improperPaths(aviews)
  J <- seq_along(aln_paths)[!file.exists(aln_paths)]
  for(j in J){
    getImproperAlignmentPairs(aviews[, j], param,
                              mapq_thr=mapq_thr,
                              use.mcols=use.mcols)
  }
  if(!all(file.exists(improperPaths(aviews)))){
    stop("Improper read pairs were not written to disk for one or more samples")
  }
  invisible()
}


thinReadPairQuery <- function(g, zoom.out=1){
  greater_20kb <- width(g) > 20e3
  g2 <- expandGRanges(g, zoom.out*width(g))
  ## for the big intervals, focus on the border areas
  if(!any(greater_20kb)){
    return(g2)
  }
  g3 <- g2[!greater_20kb]
  index <- greater_20kb
  g4 <- GRanges(seqnames(g)[index],
                IRanges(start(g)[greater_20kb]-5e3,
                        start(g)[greater_20kb]+5e3))
  g5 <- GRanges(seqnames(g)[index],
                IRanges(end(g)[greater_20kb]-5e3,
                        end(g)[greater_20kb]+5e3))
  gg <- c(granges(g3), g4, g5)
  dj <- disjoin(gg)
  dj
}


## will retrieve improper and proper reads
.mappedReadPairFlags <- function() scanBamFlag(isDuplicate=FALSE,
                                               isPaired=TRUE,
                                               isUnmappedQuery=FALSE,
                                               hasUnmappedMate=FALSE,
                                               isProperPair=TRUE)

.scan_all_readpairs <- function(granges, bam.file, flags){
  g <- reduce(granges)
  param <- ScanBamParam(flag=flags, what=c("flag", "mrnm", "mpos", "mapq"), which=g)
  ##x <- readGAlignmentsFromBam(bam.file, param=param, use.names=TRUE)
  x <- readGAlignments(bam.file, param=param, use.names=TRUE)
  x <- makeGAlignmentPairs2(x, use.mcols="flag")
  x
}

R1isFirst <- function(galp) start(first(galp)) < end(last(galp))
validFirstR1 <- function(galp)  strand(first(galp)) == "+" & strand(last(galp)) == "-" & R1isFirst(galp)
validLastR1 <- function(galp)  strand(first(galp)) == "-" & strand(last(galp)) == "+" & !R1isFirst(galp)

validPairForDeletion <- function(galp){
  ## for deletions:
  ## R1 +, R2-, R1<R2
  ## R1 -, R2+, R1>R2
  same_chrom <- chromosome(first(galp)) == chromosome(last(galp))
  ##if(!all(same_chrom)) galp <- galp[same_chrom]
  r1pos_r2neg <- validFirstR1(galp)
  r1neg_r2pos <- validLastR1(galp)
  r1pos_r2neg | r1neg_r2pos & same_chrom
}

#' Import properly paired reads from a bam file
#'
#' 
#' @return a \code{GAlignmentPairs} object
#' @export
#' @param bam_path character string providing complete path to bam file
#' @param gr a \code{GRanges} object
#' @param param a \code{DeletionParam} object
properReadPairs <- function(bam_path, gr, param=DeletionParam()){
  query <- thinReadPairQuery(gr)
  if(length(query) == 0) {
    return(.GAlignmentPairs())
  }
  seqlevelsStyle(query) <- bamSeqLevelsStyle(param)
  flags <- .mappedReadPairFlags()
  galp <- .scan_all_readpairs(query,
                              bam.file=bam_path,
                              flags=flags)
  is_valid <- validPairForDeletion(galp)
  galp <- galp[is_valid]
  galp
}

#' Extract all mapped read pairs near a deletion(s)
#'
#' Wrapper for \code{readGAlignments} that extracts mapped read pairs
#' near deletions.
#'
#' @export
#' @return a \code{GAlignmentPairs} object
#' @param object a \code{GRanges} object
#' @param bam.file a character string providing complete path to bam file
#' @param param a \code{ScanBamParam} object
readPairsNearVariant <- function(object, bam.file){
  flags <- .mappedReadPairFlags()
  ##granges <- thinReadPairQuery(object, thin)
  galp <- .scan_all_readpairs(object, bam.file=bam.file, flags=flags)
  seqlevelsStyle(galp) <- seqlevelsStyle(object)
  is_valid <- validPairForDeletion(galp)
  galp[is_valid]
}
