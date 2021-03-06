% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alignment-utils.R
\name{getImproperAlignmentPairs}
\alias{getImproperAlignmentPairs}
\title{Extract all improperly paired reads from a bam file as a GAlignmentPairs object}
\usage{
getImproperAlignmentPairs(bam.file, param = improperAlignmentParams(), build)
}
\arguments{
\item{bam.file}{complete path to BAM file}

\item{param}{a \code{ScanBamParam} object.}

\item{build}{the reference genome buld that reads were aligned to.  Currently
supported builds include "hg19" and "hg18".}
}
\value{
a \code{GAlignmentPairs} object
}
\description{
This function is a wrapper for readGAlignments. At the time this function was created, one could not create a GAlignmentPairs object of improperly paired reads using existing infrastructure in the GenomicAlignments package.  This function may be deprecated in the future.
}
\examples{
  library(svbams)
  path <- system.file("extdata", package="svbams")
  bam.file <- file.path(path, "cgov10t.bam")
  irp <- getImproperAlignmentPairs(bam.file, build="hg19")

}
\seealso{
See \code{\link[GenomicAlignments]{makeGAlignmentPairs}}
  for details regarding \code{use.mcols} argument.  See
  \code{\link{improperAlignmentParams}} for creating a
  \code{ScanBamParam} object with the appropriate flags for
  extracting improper read pairs.
}
\keyword{internal}
