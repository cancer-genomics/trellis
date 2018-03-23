#' @include help.R
NULL

#' Counts reads mapped to genomic intervals
#'
#' Utilizes the \code{countBam()} function in the \code{Rsamtools} package to
#' count number of reads mapped to ranges determined by
#' a \code{ScanBamParams()} object.
#'
#' @param object A \code{BamViews} object
#' @param scan_param An optional \code{ScanBamParam} object specifying any
#' parameters to consider when importing reads from bam files included in the
#' \code{BamViews} object.  If omitted, a \code{scanBamParams} object
#' will automatically be created to include all ranges from
#' the \code{BamViews} object and to ignore unmapped and duplicate reads.
#' @return An integer-vector of counts corresponding to each interval in the \code{BamViews} object.
#'
#' @examples
#' library(Rsamtools)
#' library(svbams)
#' library(svfilters.hg19)
#' data(bins1kb)
#' extdir <- system.file("extdata", package="svbams", mustWork=TRUE)
#' bamfile <- file.path(extdir, "cgov10t.bam")
#' bins <- keepSeqlevels(bins1kb, "chr3", pruning.mode = "coarse")
#' bins <- subsetByOverlaps(bins, GRanges("chr3", IRanges(59600000, 61000000)))
#' bviews <- BamViews(bamRanges=bins, bamPaths=bamfile)
#' counts <- binCounts(bviews)
#' head(counts)
#' length(bamRanges(bviews)) == length(counts) #TRUE
#'
#' @export
binCounts <- function(object, scan_param){
  if(missing(scan_param)){
    scan_param <- countParam(gr=bamRanges(object))
  }
  bamfile <- bamPaths(object)
  cnt <- countBam(bamfile, param=scan_param)$records
  cnt
}


#' Count reads in bamRanges of a BamViews object
#'
#' This function counts the number of reads mapped to each interval
#' given by \code{bamRanges(object)}.
#'
#' @param object A \code{BamViews} object
#' @param scan_param a \code{ScanBamParams} object
#' @export
countReads2 <- binCounts


#' Generate parameters that specify which reads in a bam file to import
#'
#' A wrapper function for \code{Rsamtools::ScanBamParam()}.
#'
#' @return A \code{\link[Rsamtools]{ScanBamParam}} object
#' @export
#' @param gr A \code{GRanges} object
#' @param flags An integer vector as returned by \code{Rsamtools::scanBamFlag()}.
#' Unmapped and duplicate reads are omitted by default.
#'
#' @examples
#' library(Rsamtools)
#' library(svfilters.hg19)
#' data(bins1kb)
#' scan_param <- countParam(bins1kb)
#' scan_param
#' class(scan_param) #ScanBamParam
countParam <- function(gr, flags=scanBamFlag(isUnmappedQuery=FALSE,
                                             isDuplicate=FALSE)){
  scan_param=ScanBamParam(flag=flags, which=gr)
}




setMethod(modev, "integer", function(x, na.rm=TRUE){
  tab <- table(x)
  as.integer(names(tab)[which.max(tab)])
})

setMethod(modev, "numeric", function(x, na.rm=TRUE){
  if(na.rm){
    notna <- !is.na(x)
    xx <- x[notna]
    dens <- density(xx)
    x_mode <- dens$x[which.max(dens$y)]
  } else {
    dens <- density(x)
    x_mode <- dens$x[which.max(dens$y)]
  }
  return(x_mode)
})

transformCounts <- function(cnt, i, centerby=c("mode", "median")){
  lx <- log2(cnt+1)
  centerby <- match.arg(centerby)
  if(centerby == "mode"){
    med <- modev(lx[i])
  } else {
    med <- median(lx[i])
  }
  lx - med
}


#' Normalize counts
#' 
#' Log-transforms the counts for each bin and centers them by by the 
#' mode (default) or median of the counts for the autosomes. 
#' @param bins a \code{GRanges} object
#' @param centerby mode or median
#' @return A numeric vector of normalized counts.
#' @details 
#' The \code{bins} argument to \code{binNormalize()} must be a \code{GRanges} object with  
#' chromosome names in the format [chr1, chr2, ...] and not [1, 2, ...] or any other naming convention.  
#' There must be a metadata column of integers named \code{cnt} corresponding to the raw counts 
#' for each bin.    
#' 
#' The output vector from \code{binNormalize()} is generally appended back to \code{bins}
#' as a metadata column named \code{std_cnt} for GC normalization by \code{binGCCorrect()}.
#' @examples
#' library(Rsamtools)
#' library(svbams)
#' library(svfilters.hg19)
#' data(bins1kb)
#' extdir <- system.file("extdata", package="svbams", mustWork=TRUE)
#' bamfile <- file.path(extdir, "cgov10t.bam")
#' bins <- keepSeqlevels(bins1kb, "chr3", pruning.mode = "coarse")
#' bins <- subsetByOverlaps(bins, GRanges("chr3", IRanges(59600000, 61000000)))
#' bviews <- BamViews(bamRanges=bins, bamPaths=bamfile)
#' bins$cnt <- binCounts(bviews)
#' head(binNormalize(bins))
#' std_cnt <- binNormalize(bins)
#' bins$std_cnt <- std_cnt
#' @export
binNormalize <- function(bins, centerby=c("mode", "median")){
  is_autosome <- chromosome(bins) %in% paste0("chr", 1:22)
  transformCounts(bins$cnt, is_autosome, match.arg(centerby))
}


isOutlier <- function(x, nmad=2) {
  mn <- median(x, na.rm=TRUE)
  sd <- mad(x, na.rm=TRUE)
  (x < mn - nmad*sd) | (x > mn + nmad*sd)
}

#' Threshold values of a numeric vector
#'
#' Threshold the values of a numeric vector according to
#' user-specified limits.  Values exceeding the threshold are reset to
#' the threshold and jittered by an amount specified by the user to
#' reduce overplotting.
#'
#' @param x numeric vector
#' @param lim numeric vector of length 2 indicating the values at which to threshold \code{x}
#' @param amount scalar indicating how much to jitter the points at the threshold.
#' @seealso \code{\link{jitter}}
#' @examples
#' x <- rnorm(10)
#' threshold(x, c(-0.5,0.5))
#' @export
threshold <- function(x, lim=c(-Inf,Inf), amount=0){
  notna <- !is.na(x)
  index1 <- which(x <= lim[1] & notna)
  index2 <- which(x >= lim[2] & notna)
  x1 <- runif(length(index1), lim[1], lim[1]+amount)
  x2 <- runif(length(index2), lim[2]-amount, lim[2])
  x[index1] <- x1
  x[index2] <- x2
  return(x)
}

fit_models_sequentially <- function(full_data, training_index, models){
  i <- training_index
  data <- full_data[i, ]
  for(j in seq_along(models)){
    fit <- eval(models[[j]])
    xhat <- predict(fit, newdata=full_data)
    ## update the full data.frame of normalized counts
    full_data$x <- full_data$x-xhat
    ## update the testing data.frame of normalized counts
    data$x <- full_data$x[i]
  }
  full_data$x
}

modelGC2 <- function(object, gc_model=.gc_model()){
  x <- assays(object)[, 1]
  is_outlier <- isOutlier(x, 2)
  is_na <- is.na(x)
  i <- which(!is_outlier & !is_na)
  if(length(i) > 50e3){
    training_index <- sample(i, min(50e3, length(i)))
  } else training_index <- i
  ##gc <- queryRanges(object)$gc
  gc <- rowRanges(object)$gc
  ## otherwise, correction will provide NAs for GC outside the range
  gc <- threshold(gc, range(gc[training_index], na.rm=TRUE))
  mp <- rowRanges(object)$map
  w <- width(rowRanges(object))
  df <- data.frame(x=x,
                   gc=gc,
                   width=w,
                   map=mp)
  x <- fit_models_sequentially(full_data=df,
                               training_index=training_index,
                               models=gc_model)
  x
}

.gc_model <- function() list(expression(loess(x~gc, data, span=1/3)))


#' Correct coverage for GC content 
#' 
#' Normalizes the coverage for a set of bins such that there is no correlation between GC content and coverage.
#' 
#' @param bins a \code{GRanges} object. This must contain a metadata column of numerics 
#' named \code{std_cnt} that is generally generated by \code{binNormalize()}. 
#' @details 
#' The relationship of counts to GC content is determined by fitting a loess curve. 
#' The lowess curve is subtracted from the counts for each bin to remove any GC bias.
#' @return A numeric vector of normalized counts for each bin.
#' @examples
#' library(Rsamtools)
#' library(svbams)
#' library(svfilters.hg19)
#' library(trellis)
#' data(bins1kb)
#' extdir <- system.file("extdata", package="svbams", mustWork=TRUE)
#' bamfile <- file.path(extdir, "cgov10t.bam")
#' bins <- keepSeqlevels(bins1kb, "chr3", pruning.mode = "coarse")
#' bins <- subsetByOverlaps(bins, GRanges("chr3", IRanges(59600000, 61000000)))
#' bviews <- BamViews(bamRanges=bins, bamPaths=bamfile)
#' bins$cnt <- trellis::binCounts(bviews)
#' std_cnt <- binNormalize(bins)
#' bins$std_cnt <- std_cnt
#' head(trellis::binGCCorrect(bins))
#' @export
binGCCorrect <- function(bins){
  std_count <- bins$std_cnt
  is_outlier <- isOutlier(std_count, nmad=2)
  is_na <- is.na(std_count)
  i <- which(!is_outlier & !is_na)
  if(length(i) > 50e3){
    training_index <- sample(i, min(50e3, length(i)))
  } else training_index <- i
  gc <- bins$gc
  ## otherwise, correction will provide NAs for GC outside the range
  gc <- threshold(gc, range(gc[training_index], na.rm=TRUE))
  mp <- bins$map
  full_data <- data.frame(x=std_count, gc=gc)
  fit <- loess(x~gc, full_data[training_index, ], span=1/3)
  xhat <- predict(fit, newdata=full_data)
  gc.adj <- full_data$x-xhat
  gc.adj
}


#' Computes a median normalized coverage across samples for each bin.
#'
#' This funciton is helpful for reducing bin-to-bin variation in normalized
#' coverage that is often correlated between samples.
#' 
#' @param nchunks The matrix of normalized coverage is potentially a very large matrix (bins x number samples).  To reduce the required RAM, we can read subsets of this matrix.  nchunks is an integer specifying how many subsets of the matrix are derived. Increasing the value of this parameter reduces the required RAM at the expense of increased computational time.
#' @param files character string of bamfile paths
#' @examples
#' library(Rsamtools)
#' library(svbams)
#' library(svfilters.hg19)
#' data(bins, package="svbams")
#' bins <- head(bins, 100)
#' extdir <- system.file("extdata", package="svbams", mustWork=TRUE)
#' bamfile <- file.path(extdir, "cgov10t.bam")
#' ## Let's pretend we have 20 BAM files
#' bamfiles <- rep(bamfile, 20)
#' tempfiles <- replicate(length(bamfiles), tempfile())
#' for(i in seq_along(bamfiles)){
#'   bviews <- BamViews(bamRanges=bins, bamPaths=bamfiles[i])
#'   bins$cnt <- binCounts(bviews)
#'   std_cnt <- binNormalize(bins)
#'   bins$std_cnt <- std_cnt
#'   gc.adj <- binGCCorrect(bins)
#'   gc.adj.int <- as.integer(round(gc.adj*1000, 0))
#'   saveRDS(gc.adj.int, file=tempfiles[i])
#' }
#' binMedians(tempfiles, nchunks=1)
#' @export
binMedians <- function(files, nchunks=50){
  ##
  ## Files are too big to load all data in memory
  ##
  ##  Autosomes + X (all samples are female): 
  ##  
  ## 1. read in a chunk for all files
  ## 2. compute medians for the chunk
  ## 3. read in all data for one sample.  Median center
  ## y_jk = x_jk - mu_j
  ##
  ##  ChrY
  ## y_jk = (x_jk - mu_j) + mu_k
  ## x = normalized coverage
  ## j = bin index
  ## k = sample index
  ## mu_j = median (across samples)
  ## mu_k = median (across bins)
  gc.adj <- readRDS(files[1])
  medians <- rep(NA, length(gc.adj))
  chunks <- sort(rep(1:nchunks, length.out=length(medians)))
  index.list <- split(seq_along(gc.adj), chunks)
  for(i in seq_along(index.list)){
    index <- index.list[[i]]
    dat <- matrix(NA, length(index), length(files))
    for(k in seq_along(files)){
      gc.adj <- readRDS(files[k])[index]
      dat[, k] <- gc.adj
    }
    meds <- matrixStats::rowMedians(dat, na.rm=TRUE)
    meds <- as.integer(round(meds, 0))
    medians[index] <- meds
  }
  medians
}

chrYMedians <- function(files, chry.index){
  ## compute the median normalized count for chromosome Y
  medians <- rep(NA, length(files))
  for(k in seq_along(files)){
    gc.adj <- readRDS(files[k])
    medians[k] <- median(gc.adj[chry.index], na.rm=TRUE)
  }
  medians
}

# readGcExperiment <- function(object, path){
#   ##path <- listDir("2gc_adj", tree)[["path"]]
#   path <- file.path(path, "2gc_adj")
#   files <- file.path(path, rdsId(object))
#   dat.list <- sapply(readRDS, files)
#   dat <- do.call(cbind, dat.list)
#   dat
# }

##--------------------------------------------------
##
## Model background
##
##--------------------------------------------------

myfilter <- function(x, filter, ...){
  res <- filter(x, filter,...)
  ##now fill out the NAs
  M <- length(filter)
  N <- (M - 1) / 2
  L <- length(x)
  for(i in seq_len(N)){
    w <- filter[(N-i+2):M]
    y <- x[1:(M-N+i-1)]
    res[i] <- sum(w*y)/sum(w)
    w <- rev(w)
    ii <- (L-(i-1))
    y <- x[(ii-N):L]
    res[ii] <- sum(w*y)/sum(w)
  }
  return(res)
}

fillInMissing <- function(x){
  if(!any(is.na(x))) return(x)
  x[1] <- 0L
  x[length(x)] <- max(x, na.rm=TRUE)
  tmp <- split(seq_along(x), x)
  tmp2 <- lapply(tmp, function(y) seq(min(y), max(y)))
  for(i in seq_along(tmp2)){
    ## these are the na indices
    j <- setdiff(tmp2[[i]], tmp[[i]])
    x[j] <- as.integer(names(tmp2)[i])
  }
  if(any(is.na(x))) stop("missing values persist")
  return(x)
}

chunkseq <- function(y){
  na.index <- which(is.na(y))
  any.na <- length(na.index) > 0
  chunks <- rep(NA, length(y))
  if(length(na.index) > 0){
    z <- y[-na.index]
  } else z <- y
  mu <- median(z)
  sigma <- mad(z)
  zlist <- split(z, sort(rep(1:20, length.out=length(z))))
  mads <- sapply(zlist, mad)
  mads2 <- median(mads)
  sigma <- min(sigma, mads2)
  condition <- (z <= (mu+2*sigma)) & (z >= (mu-2*sigma))
  k <- max(round(0.01 * length(y), 0), 5)
  smooth.condition <- myfilter(condition, rep(1/k, k))##big odd number to avoid na's
  d <- smooth.condition >= 0.5
  chks <- c(0, cumsum(abs(diff(d))))
  if(length(unique(chks)) > 4){
    chks <- rep(0, length(z))
  }
  if(any.na) {
    chunks[-na.index] <- chks
  } else chunks <- chks
  chunks <- fillInMissing(chunks)
  chunks
}

waveFun <- function(cn, span, ...){
  if(length(cn) < 100){
    bg <- rep(0, length(cn))
    return(bg)
  }
  mu <- median(cn, na.rm=TRUE)
  sigma <- mad(cn, na.rm=TRUE)
  i <- which(abs(cn - mu) <= 2*sigma)
  df <- data.frame(cn)
  df$x <- seq_len(nrow(df))
  if(nrow(df) > 10e3){
    ii <- intersect(i, seq(1, nrow(df), 2))
    fit <- loess(cn~x, df[ii, ], span=span, ...)
  } else {
    fit <- loess(cn~x, df[i, ], span=span, ...)
  }
  bg <- predict(fit, newdata=df)
  bg[is.na(bg)] <- mu
  bg
}

backgroundChunk <- function(y){
  if(length(y) > 1000) span <- 1/30 else span <- 1/5
  waveFun(y, span=span)
}

backgroundSample <- function(y){
  chunks <- chunkseq(y)
  ylist <- split(y, chunks)
  bg_list <- vector("list", length(ylist))
  for(i in seq_along(ylist)){
    bg_list[[i]] <- backgroundChunk(ylist[[i]])
  }
  bg <- do.call(c, bg_list)
  bg
}





segmentInflections <- function(object, bg){
  rowid <- rownames(object)
  names(bg) <- rowid
  dat <- CNA(bg,
             chrom=chromosome(object),
             maploc=start(object))
  seg <- segment(dat, alpha=1e-3,
                 undo.splits="sdundo",
                 verbose=0)
  gr <- cbs2granges(seg, seqinfo(object))
  gr$sample <- colnames(object)
  gr
}

inflectionFun <- function(bg){
  inflection <- rep(FALSE, length(bg))
  d <- diff(bg)
  regions <- cumsum(diff(sign(d))!=0)
  uregions <- unique(regions)
  index <- match(uregions, regions) + 2
  inflection[index] <- TRUE
  inflection[c(1, length(inflection))] <- TRUE
  inflection
}

backgroundChr <- function(se, ...){
  se <- se[, 1]
  CHR <- unique(chromosome(se))
  cn <- as.numeric(copynumber(se))
  ## fit loess through data
  bg <- backgroundSample(cn)
  inflections <- inflectionFun(bg)
  gr <- segmentInflections(se[inflections, ], bg[inflections])
  any_overlap <- overlapsAny(se, gr)
  se2 <- se[any_overlap, ]
  bg <- bg[any_overlap]
  means_list <- vector("list", length(gr))
  for(i in seq_along(gr)){
    index <- subjectHits(findOverlaps(gr[i], se2))
    means_list[[i]] <- rep(gr$seg.mean[i], length(index))
  }
  seg_means <- do.call(c, means_list)
  bg <- bg-seg_means
  result <- rep(0, length(se))
  index <- subjectHits(findOverlaps(se, se2))
  result[index] <- bg
  result
}

robustBackground <- function(se, verbose=0, ...){
  uchrom <- seqlevels(se)
  bg_list <- vector("list", length(uchrom))
  for(i in seq_along(uchrom)){
    CHR <- uchrom[i]
    se_chr <- keepSeqlevels(se, CHR, pruning.mode = "coarse")
    if(length(se_chr) < 1000){
      bg_list[[i]] <- rep(0, nrow(se_chr))
      next()
    }
    bg_list[[i]] <- backgroundChr(se_chr,
                                  verbose=verbose)
  }
  bg <- unlist(bg_list)
  bg
}

backgroundCorrect <- function(se){
  bg <- robustBackground(se)
  copynumber(se)[, 1] <- copynumber(se)[, 1] - bg
  ## apply a second time
  bg <- robustBackground(se)
  copynumber(se)[, 1] <- copynumber(se)[, 1] - bg
  se
}

modelBackground <- function(object){
  se <- as(object, "RangedSummarizedExperiment")
  tmp <- backgroundCorrect(se)
  cn <- copynumber(tmp)
  cn2 <- as.numeric(cn)
}

## DELETE
# backgroundExperiment2 <- function(object, path){
#   ##path <- listDir("3background_adj", tree)[["path"]]
#   path <- file.path(path, "3background_adj")
#   files <- file.path(path, rdsId(object))
#   J <- seq_len(ncol(object))
#   for(j in J){
#     ##xx <- rep(NA, length(ix))
#     obj <- object[, j]
#     file <- files[j]
#     x <- modelBackground(obj)
#     x <- setNames(as.integer(round(x * 1000, 0)), names(x))
#     saveRDS(x, file=file)
#   }
#   paths(object) <- files
#   setScale(object) <- 1000
#   object
# }

# readBackgroundExperiment <- function(object, path){
#   path <- file.path(path, "3background_adj")
#   files <- file.path(path, rdsId(object))
#   J <- seq_len(ncol(object))
#   stopifnot(all(file.exists(files)))
#   for(j in J){
#     ##xx <- rep(NA, length(ix))
#     obj <- object[, j]
#     file <- files[j]
#     x <- modelBackground(obj)
#     x <- setNames(as.integer(round(x * 1000, 0)), names(x))
#     saveRDS(x, file=file)
#   }
#   paths(object) <- files
#   setScale(object) <- 1000
#   object
# }


preprocessDirs <- function(){
  c("0counts", "1transformed_centered",
    "2gc_adj", "3background_adj")
}

#' Determines which rows from a \code{GAlignmentPairs} object are duplicates of a row with a smaller subscript
#'
#' Given a \code{GAlignmentPairs} object, \code{duplicatedGAlignmentPairs} returns a logical vector with one element for each row.
#'
#' @param galp a \code{GAlignmentPairs} object
#' @return a logical vector corresponding to each row in the \code{GAlignmentPairs} object
#'
#' @examples
#' library(Rsamtools)
#' library(svbams)
#' library(GenomicAlignments)
#' bam_path <- system.file("extdata", package="svbams", mustWork=TRUE)
#' bamfile <- file.path(bam_path, "cgov10t.bam")
#' data(bins, package='svbams')
#' param1 <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,
#'                                         isSecondaryAlignment=FALSE),
#'                        which=bins)
#' param1b <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,
#'                                         isSecondaryAlignment=FALSE))
#' galp1 <- readGAlignmentPairs(bamfile, param=param1)
#' galp1 <- galp1[do.call(order, as.data.frame(galp1))]
#' galp1b <- readGAlignmentPairs(bamfile, param=param1b)
#' galp1b <- subsetByOverlaps(galp1b, bins, ignore.strand=TRUE)
#' galp1b <- galp1b[do.call(order, as.data.frame(galp1b))]
#' identical(galp1, galp1b) ## FALSE because galp1 has duplicates
#' galp1 <- galp1[!duplicatedGAlignmentPairs(galp1)]
#' identical(galp1, galp1b)
#'
#' @export
duplicatedGAlignmentPairs <- function(galp){
  df <- as(galp, "data.frame")
  duplicated(df)
}

#' Counts the number of read fragments mapped to bins
#' 
#' For each bin provided by a \code{GRanges object}, \code{binFragments} 
#' counts the number of fragments, determined by read pairs in a 
#' \code{GAlignmentPair} object, that overlap that bin
#' 
#' @param galp a \code{GAlignmentPairs} object
#' @param bins a \code{GRanges} object
#' @return an integer-vector of counts corresponding to each range in the \code{GRanges} object
#' 
#' @examples 
#' library(Rsamtools)
#' library(svbams)
#' library(GenomicAlignments)
#' bam_path <- system.file("extdata", package="svbams", mustWork=TRUE)
#' bamfile <- file.path(bam_path, "cgov10t.bam")
#' data(bins, package='svbams')
#' param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,
#'                                        isSecondaryAlignment=FALSE)) 
#' galp <- readGAlignmentPairs(bamfile, param=param)
#' fragmentCounts <- binFragments(galp, bins) 
#'
#' @export
binFragments <- function(galp, bins){
  galp <- galp[!duplicatedGAlignmentPairs(galp)]
  fragments <- as(galp, "GRanges") 
  fragmentCounts <- countOverlaps(bins, fragments)
  fragmentCounts
}


