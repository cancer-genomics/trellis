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
#' topics(all_paths)
#'
#' ## Subsetting by character-string uses grep
#' all_paths[["count"]]
#' all_paths["count"]
#' all_paths["data"] ## many paths have 'data' as part of string
#' all_paths["gc"]
#'
#' @seealso \code{\link{topics}}
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

#' Lists folders for a given topic
#'
#' Lists folders for major topics, together with the path to the major
#' topic.
#'
#' TODO: add topics()
#' 
#' @examples
#' dp <- DataPaths(tempdir(), "Test")
#' listDir("preprocess", dp)
#'
#' @seealso See \code{\linkS4class{DataPaths}}.  See
#'   \code{\link{topics}} for a listing of the major headings.
#' 
#'@export
#' @param folder length-one character string
#' @param object an object of class \code{DataPaths}
listDir <- function(folder, object){
  folders <- list.files(object[folder])
  path <- object[folder]
  list(path=path, folders=folders)
}

preprocessDir <- function(object){
  file.path(dataDir(object), "preprocess")
}

segmentDir <- function(object){
  file.path(dataDir(object), "segment")
}

alignDir <- function(object){
  file.path(dataDir(object), "alignments")
}

rearDir <- function(object){
  file.path(dataDir(object), "rearrangements")
}

figDir <- function(object){
  file.path(dirname(object[1]), "figures")
}

setMethod("show", "DataPaths", function(object){
  cat("An object of class 'DataPaths'\n")
  cat(" Top level directory:\n")
  top <- dirname(object[1])
  cat(" ", top, "\n")
  cat("\n")
  cat("    /data/preprocess: \n")
  preprocess_dirs <- list.files(preprocessDir(object))
  ##preprocess_dirs <- listDir("preprocess", object)[["folders"]]
  preprocess_dirs <- paste(preprocess_dirs, collapse=" | ")  
  cat("      ", preprocess_dirs, "\n")
  ## CNVs
  cat("    /data/segments: \n")  
  cnv_dirs <- list.files(segmentDir(object))
  cnv_dirs <- paste(cnv_dirs, collapse=" | ")  
  cat("      ", cnv_dirs, "\n")
  ## alignments
  cat("    /data/alignments: \n")  
  aln_dirs <- list.files(alignDir(object))
  aln_dirs <- paste(aln_dirs, collapse=" | ")  
  cat("      ", aln_dirs, "\n")
  ## rearrangements
  cat("    /data/rearrangements: \n")  
  rear_dirs <- list.files(rearDir(object))
  rear_dirs <- paste(rear_dirs, collapse=" | ")  
  cat("      ", rear_dirs, "\n")
  ## figures
  cat("    /figures: \n")  
  fig_dirs <- list.files(figDir(object))
  fig_dirs <- paste(fig_dirs, collapse=" | ")  
  cat("      ", fig_dirs, "\n")        
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
#' topics(dp)
#' 
#' @param object a \code{DataPaths} object
#' @export
topics <- function(object){
  list.files(dataDir(object))
}

dataDir <- function(object) object[1]

dirCreate <- function(x){
  xx <- x[!file.exists(x)]
  if(length(xx) == 0) return(x)
  suppressWarnings(sapply(xx, dir.create, recursive=TRUE))
  x
}

## top-level folders
create_top_dirs <- function(path, rootname, dryrun=TRUE){
  topic_nms <- c("data", "figures",  "unit_test", "fasta", "fasta_unmapped")
  paths <- file.path(path, topic_nms)
  if(!dryrun){
    dirCreate(paths)
  }
  paths
}

create_preprocess_dirs <- function(path, rootname, dryrun=TRUE){
  topic_nms <- c("0counts", "1transformed_centered", "2gc_adj",
                 "3background_adj", "final_preprocess")
  paths <- file.path(path, file.path("data", file.path("preprocess", topic_nms)))
  if(!dryrun){
    dirCreate(paths)
  }
  paths
}

create_cnv_dirs <- function(path, rootname, dryrun=TRUE){
  topic_nms <- c("0cbs", "1deletions", "2amplicons", "final_segment")
  paths <- file.path(path, file.path("data", file.path("segment", topic_nms)))
  if(!dryrun){
    dirCreate(paths)
  }
  paths  
}

create_alignments <- function(path, rootname, dryrun=TRUE){
  topic_nms <- c("0improper", "1blat_mapped", "2blat_unmapped", "3reads", "4parsed_mapped",
                 "5parsed_unmapped")
  paths <- file.path(path, file.path("data", file.path("alignments", topic_nms)))
  if(!dryrun){
    dirCreate(paths)
  }
  paths    
}

create_rearrangements <- function(path, rootname, dryrun=TRUE){
  topic_nms <- c("0linked", "1somatic", "2blat_mapped", "3blat_unmapped", "final_rearrangement")
  paths <- file.path(path, file.path("data", file.path("rearrangements", topic_nms)))
  if(!dryrun){
    dirCreate(paths)
  }
  paths      
}


create_figures <- function(path, rootname, dryrun=TRUE){
  topic_nms <- c("0preprocess", "1dels", "2amps", "3rear", "intrachrom", "interchrom")
  paths <- file.path(path, "figures", topic_nms)
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
  rear_dirs <- create_rearrangements(path, rootname, dryrun)
  fig_dirs <- create_figures(path, rootname, dryrun)
  c(top_dirs, preprocess_dirs, cnv_dirs, aln_dirs, rear_dirs,
    fig_dirs)
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
  i <- grep(i, x)
  if(length(i) > 1){
    ## return the parent directory
    i <- i[1]
    x <- dirname(x[[i]])
  } else {
    x <- x[[i]]
  }
  x
})

#' @param j ignored
#' @param ... ignored
#' @param drop ignored
#' @rdname DataPaths-class
setMethod("[", c("DataPaths", "character"), function(x, i, j, ..., drop=FALSE){
  i <- grep(i, x)
  if(length(i) > 1){
    ## return the parent directory
    i <- i[1]
    x <- dirname(x[[i]])
  } else {
    x <- x[[i]]
  }  
  x
})
