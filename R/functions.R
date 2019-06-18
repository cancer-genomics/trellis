getBamView <- function(){
  seq.info <- seqinfo(bins1kb)
  extdata <- system.file("extdata", package="svbams")
  bview <- BamViews(bamPaths=file.path(extdata, "cgov44t_revised.bam"))
  seqinfo(bamRanges(bview)) <- seq.info
  bins <- keepSeqlevels(bins1kb, paste0("chr", c(1:22, "X")),
                        pruning.mode="coarse")
  bamRanges(bview) <- bins
  bview
}

improper.path <- function(){
  extdata <- system.file("extdata", package="svbams")
  file.path(extdata, "improper_cgov44t.bam.rds")
}

cgov44t_preprocess<- function(){
  extdata <- system.file("extdata", package="svbams")
  id <- "cgov44t_revised.bam"
  bamfile <- file.path(extdata, id)

  cnvpath <- system.file("extdata", package="svbams")
  gr <- readRDS(file.path(cnvpath, "cgov44t_segments.rds"))
  segs <- keepSeqlevels(gr, "chr15", pruning.mode="coarse")

  irp.file <- file.path(extdata, "cgov44t_improper.rds")
  irp <- readRDS(irp.file)
  ddir <- system.file("extdata", package="svbams",
                      mustWork=TRUE)
  lr <- readRDS(file.path(ddir, "preprocessed_coverage.rds"))/1000
  seqlevels(bins1kb, pruning.mode="coarse") <- paste0("chr", c(1:22, "X"))
  bins1kb$log_ratio <- lr
  del.gr <- segs[segs$seg.mean < hemizygousThr(DeletionParam())]
  reduce <- IRanges::reduce
  proper.del <- properReadPairs(bamfile,
                                gr=reduce(del.gr,
                                          min.gapwidth=2000))
  rps <- list(improper=irp, proper_del=proper.del)
  pdat <- trellis::preprocessData(bam.file=bamfile,
                                  genome=genome(segs)[[1]],
                                  segments=segs,
                                  read_pairs=rps,
                                  bins=bins1kb)
}
