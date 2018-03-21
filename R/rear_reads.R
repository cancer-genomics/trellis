#' Extracts read sequences
#'
#' This function extracts the sequences of reads from a bam file and
#' saves the result as an interemediate file.  If the intermediate
#' file already exists, the read sequences are read from disk.
#'
#' @details
#'
#' This function calls \code{getReadAlignmentPairs}, the function
#'   originally called by \code{writeImproperAlignments2}.  The reason
#'   we call this function a second time is that the read sequences
#'   were not saved in the initial query.  Saving the sequences in the
#'   original query would prevent having to access the bam file a
#'   second time, though increase the size of the improper read files.
#'   We could use.  Another improvement in efficiency might come from
#'   parsing the original large bam into a smaller bam containing
#'   improper read pairs in which both mates were mapped.
#'
#' 
#' @export
#' @param dp a \code{DataPaths} object
#' @param rlist a \code{RearrangementList}
#' @param aviews an \code{AlignmentViews2} object
#' @param MAX the maximum number of read pairs to extract for a
#'   rearrangement.  If the number of read pairs supporting a
#'   rearrangement is greater than MAX, a random sample of MAX
#'   supporting read pairs is returned.
tag_sequenceExperiment <- function(dp, rlist, aviews, MAX=25L){
  .Deprecated()
  tag_list <- vector("list", length(rlist))
  names(tag_list) <- names(rlist)
  files <- file.path(dp["3reads"], paste0(names(tag_list), ".rds"))
  for(i in seq_along(tag_list)){
    if(file.exists(files[i])){
      tags <- readRDS(files[i])
      tags <- tbl_df(tags)
      tag_list[[i]] <- tags
      next()
    }
    tags <- getSequenceOfReads(rlist[[i]], aviews[, i], MAX=MAX)
    saveRDS(tags, file=files[i])
    tag_list[[i]] <- tbl_df(tags)
  }
  tag_list
}


