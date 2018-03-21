id <- "n_PGDX4500T3_WGS_Ex.bam"
id.rds <- paste0(id, ".rds")
bamDir <- "/dcl01/scharpf/data/rscharpf/Dorothy_Ovarian2/hg18/eland/standard_header"

library(trellis)
library(Rsamtools)
library(svfilters.hg18)
devtools::load_all("/dcl01/scharpf1/data/rscharpf/dorothy/ovarian/ovariantn", export_all = FALSE)
data(bviews_hg18, package="ovariantn")  

# Loading in the preprocessed coverage for 
log_ratio <- readRDS(file.path("/dcl01/scharpf/data/rscharpf/projects/dorothy/ovarian/preprocess/data2", id.rds))
log_ratio <- log_ratio / 1000
bins <- bamRanges(bviews_hg18)
bins$log_ratio <- log_ratio

# Reading in the segmentation
segs <- readRDS(file.path("/dcl01/scharpf1/data/dbruhm/dorothy/ovarian/deletions/segmentation", id.rds))

# Reading in the improper read pairs
improper_rp <- readRDS(file.path("/dcl01/scharpf1/data/rscharpf/dorothy/ovarian/ovariantn-data/data/alignments/0improper/", id.rds))

# Only keep segments below hemizygous threshold for deletions (correct this for tumor purity)
del.gr <- reduce(segs[segs$seg.mean < hemizygousThr(DeletionParam())], min.gapwidth=2000)

# Reading in proper read pairs surrounding the putative deletions
bam.file <- file.path(bamDir, id) 
proper_rp <- properReadPairs(bam.file, gr=del.gr, DeletionParam())

# Calling deletions using sv_deletions
read_pairs <- list(proper_del=proper_rp, improper=improper_rp)
pdata <- preprocessData(bam.file=bam.file,
                        genome="hg18",
                        bins=bins, segments=segs,
                        read_pairs=read_pairs)
saveRDS(pdata, file="~/pdata.rds")
untrace(trellis:::updateImproperIndex2, browser)
trace(sv_deletions, browser)

pdata <- readRDS("pdata.rds")
trace(sv_deletions, browser)
deletions <- sv_deletions(preprocess=pdata)
