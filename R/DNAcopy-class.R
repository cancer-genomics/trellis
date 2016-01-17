#' Class for storing circular binary segmentation results
#'
#' @seealso \code{\link[DNAcopy]{segment}}
#' 
#' @exportClass DNAcopy
#' @aliases DNAcopy
#' @family DNAcopy
#' @rdname DNAcopy-methods
#' @name DNAcopy-class
setOldClass("DNAcopy")

#' As("DNAcopy", "GRanges")
#'
#' Coerce \code{DNAcopy} object to \code{GRanges}
#' 
#' @name as
#' @family DNAcopy
#'
#' @param from S3 class DNAcopy
#' @return a \code{GRanges} object
#' @examples
#' library(DNAcopy)
#' ## example from DNAcopy
#' genomdat <- rnorm(500, sd=0.1) +
#'               rep(c(-0.2,0.1,1,-0.5,0.2,-0.5,0.1,-0.2),
#'                   c(137,87,17,49,29,52,87,42))
#' chrom <- rep(1:2,c(290,210))
#' maploc <- c(1:290,1:210)
#' test1 <- segment(CNA(genomdat, chrom, maploc))
#' as(test1, "GRanges")
#' @rdname DNAcopy-methods
#' @docType methods
#' @export
setAs("DNAcopy", "GRanges", function(from, to){
  gr <- GRanges(as.character(from$output$chrom),
                IRanges(from$output$loc.start, from$output$loc.end),
                seg.mean=from$output$seg.mean)
  gr
})
