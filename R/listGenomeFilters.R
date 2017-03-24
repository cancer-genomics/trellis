#' Genomic filters for the somatic amplicon analysis
#'
#' When paired normal samples are unavailable, we use germline
#' lymphoblast cell lines to identify CNVs and outliers that are
#' common in the germline and unlikely to be somatic.  Filters
#' relevant for the somatic amplicon analysis include germline CNVs,
#' centromeric regions, and assembly gaps.  We additionally defined a
#' 'paired bin filter' containing bins that appear to be linked by
#' improperly spaced reads in the lymphoblast samples.  Aberrantly
#' spaced linked bins could indicate a new junction in the lymphoblast
#' cell line or a sequencing artifact.
#'
#' REFACTORING: Separate filters from the parameter list. Create a single filter
#' object for both the amplicon and deletion analyses. * See germline_filters.
#' Do we still need this? *
#'
#' @param genome length-one character vector specifying UCSC genome build
#' @return a named list
#'
#' @examples
#' filters <- listGenomeFilters("hg19")
#'
#' @export
#' @seealso \code{germline_filters}  
listGenomeFilters <- function(genome=c("hg19", "hg18")){
  ##if(ucsc_build != "hg19") stop("Only available for build hg19")
  genome <- match.arg(genome)
  pkg <- paste0("svfilters.", genome)
  ##data(transcripts, package=pkg, envir=environment())
  ##transcripts <- get("transcripts")
  data(assembly_gaps, package=pkg, envir=environment())
  assembly_gaps <- get("assembly_gaps")
  data(gaps, package=pkg, envir=environment())
  gaps <- get("gaps")
  centromeres <- gaps[gaps$type=="centromere"]
  seqinfo(centromeres) <- seqinfo(assembly_gaps)
  data(coverage_filters, package=pkg, envir=environment())
  coverage_filters <- get("coverage_filters")
  cnv <- reduce(c(coverage_filters[["amplicon"]],
                  coverage_filters[["deletion"]]))
  out <- reduce(coverage_filters[["outlier"]])
  list(centromeres=centromeres,
       assembly_gaps=assembly_gaps,
       germline_cnv=cnv,
       outliers=out)
       ##transcripts=transcripts)
}


#' Create a reduced set of germline and sequence filters
#'
#' @details
#'
#' The reduced set is restricted to the sequence names provided in
#'   \code{seqlev} (if provided).
#'
#' @examples
#' library(svfilters.hg19)
#' filter.list <- listGenomeFilters("hg19")
#' filter.list <- filter.list[names(filter.list) != "transcripts"]
#' filters <- reduceGenomeFilters(filter.list)
#'
#' @seealso \code{\link{listGenomeFilters}}
#'
#' @return a reduced \code{GRanges} object of germline and sequence filters
#' @export
#' @param filters a list of germline CNV filters as GRanges objects
#' @param seqlev character vector of sequence names
#' @seealso See \code{\link[GenomicRanges]{inter-range-methods}} for a
#'   description of \code{reduce}.
#'
reduceGenomeFilters <- function(filters, seqlev){
  ##filters <- listGenomeFilters()
  ##filters <- filters[-match("transcripts", names(filters))]
  r <- reduce(unlist(GRangesList(lapply(filters, granges))))
  if(missing(seqlev)) return(r)
  keepSeqlevels(r, seqlev, pruning.mode="coarse")
}

genomeFilters <- function(build=c("hg19", "hg18")){
  build <- match.arg(build)
  svfilters <- paste0("svfilters.", build)
  data(germline_filters, package=svfilters,
       envir=environment())
  reduceGenomeFilters(germline_filters)
}

#' List germline rearrangement filters derived from 10 lymphoblast
#' cell lines (mixed ethnicities) and 8 normal blood sammples.
#'
#' @param genome length-one character vector indicating genome build (hg18 or hg19)
#' @return a \code{GRangesList}
#' @examples
#' listRearFilters("hg19")
#' @seealso \code{\link{reduceRearFilters}}
#' @export
listRearFilters <- function(genome=c("hg19", "hg18")){
  genome <- match.arg(genome)
  pkg <- paste0("svfilters.", genome)
  data("lymphoblast_rear", package=pkg, envir=environment())
  lymphoblast_rear <- get("lymphoblast_rear")
  lymphoblast_rear
}

#' Provides a reduced set of germline rearrangement filters derived
#' from 10 lymphoblast cell lines (mixed ethnicities) and 8 normal
#' blood sammples.
#'
#' TODO: deprecate
#'
#' @details
#'
#' The reduced set is restricted to the sequence names provided in
#'   \code{seqlev} (if provided).
#'
#' @examples
#' filters <- listRearFilters("hg19")
#' reduceRearFilters(filters, "chr1")
#'
#' @return a \code{GRanges} object
#' @seealso \code{\link{listRearFilters}}
#' @param filters a list of germline CNV filters
#' @param seqlev character vector of sequence names
#' 
#' @export
reduceRearFilters <- function(filters, seqlev){
  ##filters <- listRearFilters()
  r <- reduce(unlist(filters))
  if(missing(seqlev)) return(r)
  keepSeqlevels(r, seqlev, pruning.mode="coarse")
}
