#' A class for storing file paths to intermediate data types
#'
#' @export
setClass("DataPaths", contains="character")

#' A class for storing file paths to intermediate data types
#'
#' Intermediate data is saved to disk to avoid rerunning long-steps in
#' computation. The Paths constructor initializes a vector of file
#' paths for storing different types of intermediate data.
#'
#' @examples
#' ## by default, file directories are not created
#' DataPaths()
#' ## Create directories for intermediate data types
#' all_paths <- DataPaths(dryrun=FALSE)
#' all_paths
#' ## Just list the 'data' subdirectories
#' dataTypes(all_paths)
#'
#' @export
#' @param path A length-one character vector specifying where the top-level directory should reside
#' @param rootname A length-one character vector providing the name of the directory to create under \code{path}
#' @param dryrun logical.  If \code{FALSE}, all directories are
#'   created. Otherwise, only the character-vector is returned.
#' @rdname DataPaths-class
DataPaths <- function(path=tempdir(), rootname="sv_data", dryrun=TRUE){
  obj <- projectTree(path, rootname, dryrun)
  as(obj, "DataPaths")
}

setValidity("DataPaths", function(object){
  msg <- TRUE
  if(!file.exists(object[["top"]])){
    msg <- "top-level directory must exist"
    return(msg)
  }
})

setMethod("show", "DataPaths", function(object){
  cat("An object of class 'DataPaths'\n")
  cat("  Top-level directory:", object[["top"]], "\n")
  cat("  see dataTypes()\n")
})

#' List directories of intermediate data types
#'
#' @seealso \code{\link{DataPaths}}
#' @return a \code{character}-vector
#' @examples
#' dirs <- DataPaths()
#' basename(dataTypes(dirs))
#' 
#' @param object a \code{DataPaths} object
#' @export
dataTypes <- function(object) object[["data"]]

dirCreate <- function(x){
  xx <- x[!file.exists(x)]
  if(length(xx) == 0) return(x)
  sapply(xx, dir.create, recursive=TRUE)
  x
}

projectTree <- function(path, rootname, dryrun=TRUE){
  topic_nms <- c("top", "data", "figures", "extdata", "versions",
                 "unit_test", "segments", "counts", "normalized_counts",
                 "gc_adjust",
                 "M",
                 "filters",
                 "fig:genome",
                 "improper",
                 "rear:improper",
                 "proper",
                 "deletions",
                 "rear:candidates",
                 "rear:filter",
                 "rear:seq",
                 "fasta",
                 "fasta:unmapped",
                 "blat:unmapped",
                 "blat:parsed",
                 "fusions",
                 "fig:ambiguous",
                 "fig:unambiguous",
                 "export",
                 "blat:alignment",
                 "amp:graphs",
                 "fig:circos",
                 "fig:rear",
                 "rearrangements",
                 "rear_stats",
                 "rearranged_reads")
  dirnames <- topic_nms
  dirnames[1] <- file.path(path, rootname)
  ## subfolders
  index <- 4:length(topic_nms)
  dirnames[index] <- c("inst/extdata", "versions", "data/unit_tests",
                       "data/segments", "data/counts",
                       "data/normalized_counts",
                       "data/gc_adjust",
                       "data/M",
                       "inst/extdata/filters",
                       "figures/genome",
                       "data/improper",
                       "data/improper_rearrangement",
                       "data/proper",
                       "data/deletions",
                       "data/rear_candidates",
                       "data/rear_filtered",
                       "data/rear_sequences",
                       "data/fasta",
                       "data/fasta-unmapped",
                       "data/blat_alignment-unmapped",
                       "data/blat_parsed",
                       "data/fusions",
                       "figures/ambiguous",
                       "figures/unambiguous",
                       "data/export",
                       "data/blat_alignment",
                       "data/amplicon_graphs",
                       "figures/circos",
                       "figures/rear",
                       "data/rearrangements",
                       "data/rearrangement_stats",
                       "data/rearranged_reads")
  dirnames[-1] <- file.path(dirnames[1], dirnames[-1])
  dirnames <- setNames(dirnames, topic_nms)
  ix <- order(names(dirnames))
  dirnames <- dirnames[ix]
  if(dryrun){
    return(dirnames)
  }
  dirCreate(dirnames)
}

