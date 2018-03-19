#' Compute the autocorrelation at a specific lag
#'
#' A wrapper for \code{acf} that returns a single autocorrelation for
#' a specified lag.
#' 
#' @seealso \code{\link{acf}}
#' @export
#' @param x a numeric vector
#' @param lag.max integer specifying the autocorrelation lag
#' @param type character string
#' @param plot logical
#' @param na.action function to handle missing observations
#' @param demean logical
#' @param ... additional arguments passed to \code{acf}
acf2 <- function(x, lag.max=10, type = c("correlation", "covariance", "partial"),
                 plot = FALSE, na.action = na.omit, demean = TRUE,
                 ...){
  y <- acf(x, lag.max=lag.max, type=type, plot=plot,
             na.action=na.action, demean=demean, ...)
  y <- y[[1]][lag.max+1, , 1]
}
