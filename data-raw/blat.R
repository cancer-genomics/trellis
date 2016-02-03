library(svalignments)
library(normalblood)
dp <- projectBlood()
data(ids)
id <- ids[4]

blat_file <- file.path(dp[["2blat_unmapped"]], paste0(id, ".txt"))
blat <- readBlat(blat_file)
rlist <- readRDS(file.path(dp["final_rearrangement"], paste0(id, ".rds")))
lb <- linkedBins(rlist)
## focus on interchromosomal
is_interchrom <- chromosome(lb) != chromosome(lb$linked.to)
rlist <- rlist[is_interchrom]
rearrangement_list <- rlist
save(rearrangement_list, file="~/Software/svpackages/svalignments/data/rearrangement_list.rda")

is_na <- is.na(blat$Tstart)
if(any(is_na)){
  blat <- blat[!is_na, ]
}
blat_gr <- GRanges(blat$Tname, IRanges(blat$Tstart, blat$Tend),
                   match=blat$match, qname=blat$Qname,
                   qstart=blat$Qstart, qend=blat$Qend)
lb <- linkedBins(rlist)
##
## A blat record must overlap one of the intervals in a linked bin
##
maxgap <- 500
record_overlaps <- overlapsAny(blat_gr, lb, maxgap=maxgap) |
  overlapsAny(blat_gr, lb$linked.to, maxgap=maxgap)
records <- blat_gr[record_overlaps]

blat_unmapped <- blat[blat$Qname %in% records$qname, ]
save(blat_unmapped, file="~/Software/svpackages/svalignments/data/blat_unmapped.rda")
