#' @include help.R
NULL

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
#' ## Subsetting by character-string uses grep
#' all_paths[["count"]]
#' all_paths["count"]
#' all_paths["data"] ## many paths have 'data' as part of string
#' all_paths["gc"]
#'
#' @export
#' @param path A length-one character vector specifying where the top-level directory should reside
#' @param rootname A length-one character vector providing the name of the directory to create under \code{path}
#' @param dryrun logical.  If \code{FALSE}, all directories are
#'   created. Otherwise, only the character-vector is returned.
#' @rdname DataPaths-class
DataPaths <- function(path=tempdir(), rootname="sv_data", dryrun=FALSE){
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

#' List folders at a particular level in the directory tree
#'
#' @examples
#' dp <- DataPaths(tempdir(), "Test")
#' listDir("preprocess", dp)
#'
#' @seealso See \code{\linkS4class{DataPaths}}.  See
#'   \code{\link{dataTypes}} for a listing of the major headings.
#' 
#'@export
#' @param folder length-one character string
#' @param object an object of class \code{DataPaths}
listDir <- function(folder, object){
  ##ix <- grep(pattern, object)
  ix <- grep(folder, basename(object))[1]
  if(is.na(ix)){
    ix <- grep(folder, basename(dirname(object)))[1]
    paths <- list(path=object[ix],
                  folders=list.files(dirname(object)[ix]))
    return(paths)
  }
  if(is.na(ix)){
    stop("pattern does not match basenames")
  }
  list(path=object[ix],
       folders=list.files(object[ix]))
}

#' Accessor for data directory
#'
#' @examples
#' dp <- DataPaths(rootname="test")
#' dataDir(dp)
#' @export
#' @param object a \code{DataPaths} object
dataDir <- function(object){
  listDir("data", object)[["path"]]
}

setMethod("show", "DataPaths", function(object){
  cat("An object of class 'DataPaths'\n")
  cat(" Top level directory:\n")
  top <- dirname(object[1])
  cat(" ", top, "\n")
  cat("\n")
  cat("    /data/preprocess: \n")
  preprocess_dirs <- listDir("preprocess", object)[["folders"]]
  preprocess_dirs <- paste(preprocess_dirs, collapse=" | ")  
  cat("      ", preprocess_dirs, "\n")
  ## CNVs
  cat("    /data/segments: \n")  
  cnv_dirs <- listDir("segment", object)[["folders"]]
  cnv_dirs <- paste(cnv_dirs, collapse=" | ")  
  cat("      ", cnv_dirs, "\n")
  ## alignments
  cat("    /data/alignments: \n")  
  aln_dirs <- listDir("alignments", object)[["folders"]]
  aln_dirs <- paste(aln_dirs, collapse=" | ")  
  cat("      ", aln_dirs, "\n")    
  cat("  see ?listDir\n")
})

top <- function(object){
  dirname(object[1])
}

#' List directories of intermediate data pipelines
#'
#' @seealso \code{\link{DataPaths}}
#' @return a \code{character}-vector
#' @examples
#' dp <- DataPaths()
#' dataTypes(dp)
#' 
#' @param object a \code{DataPaths} object
#' @export
dataTypes <- function(object) {
  d <- listDir("data", object)
  d[["folders"]]
}

dirCreate <- function(x){
  xx <- x[!file.exists(x)]
  if(length(xx) == 0) return(x)
  suppressWarnings(sapply(xx, dir.create, recursive=TRUE))
  x
}

## top-level folders
create_top_dirs <- function(path, rootname, dryrun=TRUE){
  topic_nms <- c("data", "figures", "extdata", "versions",
                 "unit_test")
  paths <- file.path(path, topic_nms)
  if(!dryrun){
    dirCreate(paths)
  }
  paths
}

create_preprocess_dirs <- function(path, rootname, dryrun=TRUE){
  topic_nms <- c("0counts", "1transformed_centered", "2gc_adj",
                 "3background_adj")
  paths <- file.path(path, file.path("data", file.path("preprocess", topic_nms)))
  if(!dryrun){
    dirCreate(paths)
  }
  paths
}

create_cnv_dirs <- function(path, rootname, dryrun=TRUE){
  topic_nms <- c("0cbs", "1deletions", "2amplicons")
  paths <- file.path(path, file.path("data", file.path("segment", topic_nms)))
  if(!dryrun){
    dirCreate(paths)
  }
  paths  
}

create_alignments <- function(path, rootname, dryrun=TRUE){
  topic_nms <- c("0improper")
  paths <- file.path(path, file.path("data", file.path("alignments", topic_nms)))
  if(!dryrun){
    dirCreate(paths)
  }
  paths    
}

projectTree <- function(path, rootname, dryrun=TRUE){
  path <- file.path(path, rootname)
  top_dirs <- create_top_dirs(path, rootname, dryrun)
  preprocess_dirs <- create_preprocess_dirs(path, rootname, dryrun)
  cnv_dirs <- create_cnv_dirs(path, rootname, dryrun)
  aln_dirs <- create_alignments(path, rootname, dryrun)
  c(top_dirs, preprocess_dirs, cnv_dirs, aln_dirs)
}

unitTestDir <- function(object) {
  object[["top"]][grep("unit_test", object[["top"]])]
}

figureDir <- function(object){
  object[["top"]][grep("unit_test", object[["figures"]])]
}

#' @param x a \code{DataPaths} object
#' @param i a length-one character vector
#' @rdname DataPaths-class
setMethod("[[", c("DataPaths", "character"), function(x, i){
  i <- grep(i, x)[1]
  x[[i]]
})

#' @param j ignored
#' @param ... ignored
#' @param drop ignored
#' @rdname DataPaths-class
setMethod("[", c("DataPaths", "character"), function(x, i, j, ..., drop=FALSE){
  i <- grep(i, x)
  x[i]
})
