#' Load Txdb object
#'
#' @param build genome build (one of hg18, hg19, or hg38)
#' @return a \code{TxDb} object
#' @examples
#' library(svfilters.hg19)
#' tx <- loadTx("hg19")
#' tx
#' @export
loadTx <- function(build){
  pkg <- paste0("svfilters.", build)
  data("transcripts", envir=environment(), package=pkg)
  tx <- get("transcripts")
  tx
}
