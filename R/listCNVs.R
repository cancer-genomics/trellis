#' List amplicons and deletions identified for a single sample
#'
#' @return a list
#' @export
#' @param amplicon_graph an \code{AmpliconGraph} object
#' @param del_sv a \code{StructuralVariant} object from the deletion pipeline
listCNVs <- function(amplicon_graph, del_sv){
  ##amplicon_graphs <- readRDS(file)
  ##amplicon_graphs <- readRDS(file.path(dp["2amplicons"], paste0(id, ".rds")))
  amplicons <- ampliconRanges(amplicon_graph)
  ##deletion_sv <- readRDS(file.path(dp["1deletions"], paste0(id, ".rds")))
  deletions <- variant(del_sv)
  list(amplicons=amplicons,
       deletions=deletions)
}
