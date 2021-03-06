% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rearrangement-utils.R
\name{rearDataFrame}
\alias{rearDataFrame}
\title{Create a data.frame of rearranged read pairs and split reads supporting a rearrangement}
\usage{
rearDataFrame(r, build, maxgap = 5000)
}
\arguments{
\item{r}{a \code{Rearrangement} object}

\item{build}{character string indicating genome build (only hg19 and hg18 currently supported)}

\item{maxgap}{the allowable distance between a split or paired end read and a transcript for assessing whether coding regions}
}
\description{
This function is useful for converting a Rearrangement object to a data.frame for subsequent visualization by ggplot.
}
\examples{
  extdata <- system.file("extdata", package="svbams")
  rfile <- file.path(extdata, "CGOV11T_1.bam.rds")
  rlist <- readRDS(rfile)
  r <- rlist[[1]]
  r2 <- fiveTo3Prime(r, "hg19")
  rearDataFrame(r2[[1]], "hg19") 
}
