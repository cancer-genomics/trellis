library(svbams)
library(svfilters.hg19)
library(svcnvs)
library(svalignments)
library(svrearrange)
extdata <- system.file("extdata", package="svbams")
id <- "cgov44t_revised.bam"
bam.file <- file.path(extdata, id)
irp.file <- file.path(extdata, "cgov44t_improper.rds")
irp <- readRDS(irp.file)
seqlevels(bins1kb, pruning.mode="coarse") <- paste0("chr", c(1:22, "X"))
seqlevels(irp, pruning.mode="coarse") <- seqlevels(bins1kb)
data(amplicon_graph, package="svcnvs")
data(deletions2, package="svcnvs")
genome_filters <- reduceGenomeFilters(germline_filters, 
                                      seqlevels(bins1kb))
filters <- rFilters(amplicons=ampliconRanges(amplicon_graph),
                    deletions=variant(deletions2),
                    rear=germline_rear,
                    germline=genome_filters)
rdat <- rearrangementData(bins=bins1kb,
                          read_pairs=list(improper=irp),
                          filters=filters)
rparam <- RearrangementParams(min_number_tags_per_cluster=5,
                              rp_separation=10e3)
minNumberTagsPerCluster(rparam)
rpSeparation(rparam)
rlist <- findCandidates2(rdat, rparam)
rdat$rlist <- rlist
rlist <- filterRear(rdat, rparam)

set.seed(123)
tags <- getSequenceOfReads(rlist, bam.file, MAX=25L)
dir <- tempdir()
fa.file <- file.path(dir, paste0(id, ".fa"))
writeToFasta(tags, fa.file)

blatpath <- "~/bin/x86_64/blat"
refgenome <- "~/Dropbox/reference_genome/hg19.fa"
outfile <- tempfile()
cmd <- paste(blatpath, refgenome, fa.file, outfile)
system(cmd)
file.copy(outfile, "../inst/extdata/blat_alignment.txt")


query <- uncouple(linkedBins(rlist))
unmapped <- unmapped_read(bam.file, query, yield_size=200000)
length(unmapped)
dir <- tempdir()
mapped_unmapped.fa <- file.path(dir, "mapped-unmapped.fa")
writeUnmappedToFasta(unmapped, mapped_unmapped.fa)

outfile <- tempfile()
cmd <- paste(blatpath, refgenome, mapped_unmapped.fa, outfile)
system(cmd)
file.copy(outfile, "../inst/extdata/blat_unmapped.txt")

blat <- readBlat("../inst/extdata/blat_unmapped.txt")
rearranged_reads <- rearrangedReads(rlist, blat, 500)
saveRDS(rearranged_reads, file="../inst/extdata/cgov44t_split_reads.rds")
