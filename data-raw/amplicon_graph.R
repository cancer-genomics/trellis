library(GenomicRanges)
library(Rsamtools)
library(svcnvs)
library(svfilters.hg19)
library(svalignments)

ddir <- system.file("extdata", package="svpreprocess", mustWork=TRUE)
cov.file <- file.path(ddir, "preprocessed_coverage.rds")
log_ratio <- readRDS(cov.file)/1000
data(bins1kb, package="svfilters.hg19")
seqlevels(bins1kb, pruning.mode="coarse") <- paste0("chr", c(1:22, "X"))
bins1kb$log_ratio <- log_ratio

path <- system.file("extdata", package="svcnvs")
segments <- readRDS(file.path(path, "cgov44t_segments.rds"))

extdata <- system.file("extdata", package="svbams")
bam.file <- file.path(extdata, "cgov44t_revised.bam")
iparams <- improperAlignmentParams(what=c("flag", "mrnm", "mpos", "mapq"))
improper_rp <- getImproperAlignmentPairs(bam.file,
                                         param=iparams)
segs <- segments
del.gr <- reduce(segs[segs$seg.mean < hemizygousThr(DeletionParam())],
                 min.gapwidth=2000)
proper_rp <- properReadPairs(bam.file, gr=del.gr, DeletionParam())
read_pairs <- list(proper_del=proper_rp, improper=improper_rp)

pdata <- preprocessData(bam.file=bam.file, genome="hg19",
                        bins=bins1kb,
                        segments=segments, read_pairs=read_pairs)
params <- ampliconParams()
amplicon_graph <- svcnvs:::sv_amplicons2(pdata, params=params)
save(amplicon_graph, file="../data/amplicon_graph.rda")
