#' @include Deletion-class.R
NULL


##--------------------------------------------------
##
## constructors
##
##--------------------------------------------------

#' Constructor for an empty GAlignmentPairs object
#' 
#' @return an empty \code{GAlignmentPairs} object
#' @export
#' @param first a \code{GAlignments}
#' @param last a \code{GAlignments}
#' 
#' @param isProperPair logical for whether to require that the reads
#'   be in a proper pair
.GAlignmentPairs <- function() GAlignmentPairs(first=GAlignments(),
                                               last=GAlignments(),
                                               isProperPair=logical())
