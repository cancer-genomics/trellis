% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/copynumber-methods.R
\docType{methods}
\name{copynumber}
\alias{copynumber}
\alias{copynumber<-}
\alias{copynumber,RangedSummarizedExperiment-method}
\alias{copynumber<-,RangedSummarizedExperiment,matrix-method}
\alias{copynumber,StructuralVariant-method}
\alias{copynumber<-,StructuralVariant,ANY-method}
\title{Accessor for 'copy' assays}
\usage{
copynumber(object)

copynumber(object) <- value

\S4method{copynumber}{RangedSummarizedExperiment}(object)

\S4method{copynumber}{RangedSummarizedExperiment,matrix}(object) <- value

\S4method{copynumber}{StructuralVariant}(object)

\S4method{copynumber}{StructuralVariant,ANY}(object) <- value
}
\arguments{
\item{object}{A \code{RangedSummarizedExperiment}}

\item{value}{matrix of copy number estimates}
}
\description{
Extract matrix of log2-transformed estimates of copy number
relative to autosomal mode or median
}
