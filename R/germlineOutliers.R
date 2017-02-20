outlierFrequencyPerBin <- function(pviews, NMAD=5){
  is_outlier <- matrix(NA, nrow(pviews), ncol(pviews))
  for(j in seq_len(ncol(pviews))){
    cat(".")
    cn <- assays(pviews[, j])[, 1]
    sd <- mad(cn, na.rm=TRUE)
    is_outlier[, j] <- cn > NMAD*sd | cn < -NMAD*sd
  }
  freq <- rowSums(is_outlier)
  freq
}

#' Identify genomic regions with outlier preprocess read depth
#' estimates from the lymphoblast cell lines
#'
#' For each 1kb bin along the genome, we assess whether two or more
#' ovarian samples have a preprocessed read depth estimate that is
#' more than \code{NMAD}s from zero.
#'
#' REFACTOR: Each major directory in the DataPaths should have a final
#' subdirectory with a views object.
#'
#' @export
#' @return a \code{GRanges} object of the reduced outlier genomic intervals
#' @param pviews a \code{PreprocessViews2} object
#' @param NMAD a length-one numeric vector indicating the number of
#'   mads from zero
germlineOutliers <- function(pviews, NMAD=5){
  freq <- outlierFrequencyPerBin(pviews, NMAD=NMAD)
  outlier_bins <- rowRanges(pviews)[freq >= 2]
  reduce(outlier_bins)
}
