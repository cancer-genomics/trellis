% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fusion-utils.R
\name{inFrameNoStop}
\alias{inFrameNoStop}
\title{Combine the in-frame and no-premature-stop LogicalLists}
\usage{
inFrameNoStop(nostop.list, inframe.list)
}
\arguments{
\item{nostop.list}{a length=3 \code{List} of \code{LogicalList} elements}

\item{inframe.list}{a length=3 \code{List} of \code{LogicalList} elements}
}
\description{
This function combines two \code{LogicalList}s for each reading frame, creating a new \code{LogicalList} list that is \code{TRUE} if the rearrangement is in-frame and has no premature stop codons.
}
\seealso{
\code{\link{inFrameList}} \code{\link{noPrematureStop}}
}
