###
install.packages("../trellis/", repos = NULL, type="source")
library(Rsamtools)
library(magrittr)
library(Rsamtools)
library(ggplot2)
library(svbams)
library(svfilters.hg19)
library(tidyverse)
library(gridExtra)
library(trellis)
library(ggnet)
library(graph)
library(Rgraphviz)
library(network)
library(sna)
library(kableExtra)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.refGene)
library(devtools)
library(testthat)
library(seqinr)
load_all("../trellis")
#bamfile="/dcl01/scharpf/data/rscharpf/OvarianAlignments/hg19/eland/CGOV44T.bam"
##Code from Vignette
extdata <- system.file("extdata", package="svbams")
bamfile <- file.path(extdata, "cgov44t_revised.bam")
what <- c("flag", "mrnm", "mpos", "mapq")
iparams <- improperAlignmentParams(what=what)
improper_rp <- getImproperAlignmentPairs(bamfile,
                                         param=iparams,
                                         build="hg19")

data(bins1kb, package="svfilters.hg19")
bins <- keepSeqlevels(bins1kb, "chr15",
                      pruning.mode="coarse")
flags <- scanBamFlag(isDuplicate=FALSE,
                     isPaired=TRUE,
                     isUnmappedQuery=FALSE,
                     hasUnmappedMate=FALSE,
                     isProperPair=TRUE)
bviews <- BamViews(bamPaths=bamfile, bamRanges=bins)
bins$cnt <- binnedCounts(bviews)

bins <- bins[ bins$cnt > 0 ]

bins$std_cnt <- binNormalize(bins)
set.seed(123)
bins$log_ratio <- binGCCorrect(bins)

g <- segmentBins(bins, param=SegmentParam())
g

data(bins1kb, package="svfilters.hg19")
ddir <- system.file("extdata", package="svbams",
                    mustWork=TRUE)
lr <- readRDS(file.path(ddir, "preprocessed_coverage.rds"))/1000
seqlevels(bins1kb, pruning.mode="coarse") <- paste0("chr", c(1:22, "X"))
bins1kb$log_ratio <- lr
bins <- keepSeqlevels(bins1kb, c("chr5", "chr8", "chr15"),
                      pruning.mode="coarse")
path <- system.file("extdata", package="svbams")
segs <- readRDS(file.path(path, "cgov44t_segments.rds"))
seqlevels(segs, pruning.mode="coarse") <- seqlevels(bins)

segs.df <- as_tibble(segs)

dp <- DeletionParam(remove_hemizygous=FALSE)
dp
del.gr <- IRanges::reduce(segs[segs$seg.mean < hemizygousThr(dp)],
                          min.gapwidth=2000)
proper_rp <- properReadPairs(bamfile, gr=del.gr, dp)
improper_rp <- keepSeqlevels(improper_rp, seqlevels(segs),
                             pruning.mode="coarse")
read_pairs <- list(proper_del=proper_rp, improper=improper_rp)

pdata <- preprocessData(bam.file=bamfile,
                        genome="hg19",
                        bins=bins,
                        segments=segs,
                        read_pairs=read_pairs)
deletions <- sv_deletions(preprocess=pdata, param=dp)
saveRDS(deletions, "data/sv_deletions.rds")
deletions <- readRDS("data/sv_deletions.rds")

###Code for full BAM file
bamfile <- "../../../Documents/Lab/Ovarian/Bamfiles/LP6005392-DNA_C03.bam"
query <- GRanges(proper(deletions[4]))
unmapped <- unmapped_read(bamfile, query, yield_size=200000)
length(unmapped)
dir <- tempdir()
saveRDS(unmapped, "data/cgov44t_del_unmapped.rds")
writeUnmappedToFasta(unmapped, "data/cgov44t_mapped_unmapped.fa")
write.table(unmapped$snms,"../../Lab/Ovarian/pancreatic_cancer_cell_lines-2019-09-23/data/muids.txt")
mapped_unmapped.fa <- file.path("data/cgov44t_mapped_unmapped.fa")

blatpath <- "~/bin/blat/blat"
refgenome <- "~/Documents/Lab/Ovarian/pancreatic_cancer_cell_lines-2019-09-23/data/hg19.fa"

outfile <- "data/blat_unammped.txt"
cmd <- paste(blatpath, refgenome, mapped_unmapped.fa, outfile)
system(cmd)
file.copy(outfile,"../../svalignments/inst/extdata/blat_unmapped.txt")

###on the cluster
mapped_unmapped.fa <- file.path("data/cgov44t_mapped_unmapped.fa")

blatpath <- "~/bin/blat"
refgenome <- "/dcl01/scharpf1/data/pipeline_resources/etc/hg19/genome/hg19.fa"

outfile <- "data/blat_unmapped.txt"
cmd <- paste(blatpath, refgenome, mapped_unmapped.fa, outfile)
system(cmd)
file.copy(outfile,"../../svalignments/inst/extdata/blat_unmapped.txt")

##

unmap.file <- file.path("../../Lab/Ovarian/pancreatic_cancer_cell_lines-2019-09-23/data/blat_unmapped.txt")
blat_unmap <- readBlat(unmap.file)
blat_unmap <- blat_unmap[rownames(unique(blat_unmap)),]
blat_unmap  <- blat_unmap[blat_unmap$Tname == "chr15",]
granges(splitReads(deletions[4]))
is_overlap()

linkedTo(rlist)[1]

split_reads <- rearrangedReads(linkedBins(rlist), blat_unmap, 500)
elementNROWS(split_reads)
split_reads
splitReads(rlist) <- split_reads
is_valid <- is_valid_splits(rlist, maxgap=50)
at_least_one <- lengths(splitReads(rlist)) >= 1
rlist2 <- rlist[ is_valid & at_least_one ]
grange