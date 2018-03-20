id <- "CGOV44T.bam"
id.rds <- paste0(id, ".rds")
ddir <- "~/Dropbox/OvarianCellLines"
sv <- readRDS(file.path(ddir, "structuralvar/data/segment/1deletions", id.rds))
sv <- sv["sv21"]
deletion <- sv
save(deletion, file="../data/deletion.rda")
