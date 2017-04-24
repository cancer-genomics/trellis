library(rmarkdown)
render("../../svrearrange/vignettes/svrearrange.Rmd")
saveRDS(rlist, file="../inst/extdata/rlist_cgov44t.rds")
