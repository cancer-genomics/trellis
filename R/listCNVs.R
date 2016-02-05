#' List amplicons and deletions identified for a single sample
#'
#' @return a list
#' @export
#' @param dp a \code{DataPaths} object
#' @param id length-one character vector providing sample id
listCNVs <- function(dp, id){
  amplicon_graphs <- readRDS(file.path(dp["2amplicons"], paste0(id, ".rds")))
  amplicons <- ampliconRanges(amplicon_graphs)
  deletion_sv <- readRDS(file.path(dp["1deletions"], paste0(id, ".rds")))
  deletions <- variant(deletion_sv)
  list(amplicons=amplicons,
       deletions=deletions)
}
