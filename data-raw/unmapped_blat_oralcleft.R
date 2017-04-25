fudir <- "~/Dropbox/Labs/Ruczinski/Fu2017/split_reads"
source(file.path(fudir, "scripts/mapped-unmapped.R"))
source(file.path(fudir, "scripts/blat_map-unmap.R"))
file.copy(file.path(fudir, "data/blat_mapped-unmapped.txt"),
          "../inst/extdata/blat_mapped-unmapped-oralcleft.txt")
