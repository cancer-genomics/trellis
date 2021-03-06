% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/germlineOutliers.R
\name{germlineOutliers}
\alias{germlineOutliers}
\title{Identify genomic regions with outlier preprocess read depth
estimates from the lymphoblast cell lines}
\usage{
germlineOutliers(pviews, NMAD = 5)
}
\arguments{
\item{pviews}{a \code{PreprocessViews2} object}

\item{NMAD}{a length-one numeric vector indicating the number of
mads from zero}
}
\value{
a \code{GRanges} object of the reduced outlier genomic intervals
}
\description{
For each 1kb bin along the genome, we assess whether two or more
ovarian samples have a preprocessed read depth estimate that is
more than \code{NMAD}s from zero.
}
\details{
REFACTOR: Each major directory in the DataPaths should have a final
subdirectory with a views object.
}
