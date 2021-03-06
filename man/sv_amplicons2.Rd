% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/amplicon-utils.R
\name{sv_amplicons2}
\alias{sv_amplicons2}
\title{Construct an AmpliconGraph from a BamViews object}
\usage{
sv_amplicons2(preprocess, amplicon_filters, params = ampliconParams())
}
\arguments{
\item{preprocess}{a list of the preprocessed data (see \code{preprocessData})}

\item{amplicon_filters}{a \code{GRanges} object of germline filters.  If this
argument is not specified then a default set of germline filters will be used
from the \code{svfilters.hg19} package if the \code{genome} slot in \code{preprocess} is
set to "hg19" or from the \code{svfilters.hg18} package if it's set to "hg18".}

\item{params}{a list of parameters for the amplicon analysis.  Default parameters
are given by \code{AmpliconParams()}.}
}
\value{
an \code{AmpliconGraph} object
}
\description{
This function constructs an \code{AmpliconGraph} from a
\code{BamViews} object for a single sample.  By default, the seeds
of the graph are focal amplicons with fold-change of nearly 3
relative to the diploid genome (log2(2.75)).  The threshold of
seeding amplicons can be adjusted by the \code{ampliconParams}
function. After seeding the graph with high-copy focal amplicons,
both neighboring (flanking) and distant low-copy focal amplicons
are added to the graph object.  Next, improperly paired reads in
which both the first and last read align to any queryRange of the
graph object are parsed from the bam file ( this function will
throw an error if not all files in \code{bamPaths} exist).  If 5 or
more improperly paired reads bridge a node to another node, these
amplicons are grouped.  Further, if a low-copy amplicon is bridged
to an existing node, the low-copy amplicon will become a node in
the graph. Amplicon groups are defined by the edges between nodes,
where the edges represent improperly paired reads that support a
junction between two amplicons.
}
\examples{
data(pdata)
pdata$bam.file <- system.file("extdata", "cgov44t_revised.bam",
                              package="svbams")
sv_amplicons2(pdata)
}
\seealso{
See \code{ampliconParams} for default parameters. The
  wrapper \code{\link{sv_amplicon_exp}} constructs and saves an
  \code{\linkS4class{AmpliconGraph}} for each sample in an
  experiment.
}
