#' Load Txdb object
#'
#' @return a \code{TxDb} object
#' @examples
#' library(svfilters.hg19)
#' tx <- loadTx()
#' tx
#' @export
loadTx <- function(build){
  pkg <- paste0("svfilters.", build)
  data(transcripts, envir=environment(), package=pkg)
  transcripts
}
