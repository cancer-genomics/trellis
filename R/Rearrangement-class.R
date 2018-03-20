#' @include GAlignmentPairs-utils.R
NULL

#' A class for storing rearrangement data
#'
#' This class contains the genomic intervals corresponding to improper
#' read clusters.  Improper read clusters that are linked are
#' candidate rearrangements.  The linking is represented by a
#' \code{GRanges} object with a 'linked.to' variable in the
#' \code{mcols}.  All improper reads belonging to a cluster,
#' irrespective of whether they link two clusters, are included and
#' indexed to facilitate subsetting operations.
#'
#' @slot linkedBins a \code{GRanges} object with 'linked.to' in the
#'   element metadata (\code{mcols})
#'
#' @slot modal_rearrangement character string indicating type of rearrangement
#'
#' @slot improper a \code{GAlignmentPairs} object of the improper read pairs
#'
#' @slot partitioning integer vector
#' @slot link1id character string label for first read cluster
#' @slot link2id character string label for second read cluster
#' @slot tags a \code{GRanges} object of single tags
#'
#' @slot tag_map_to_linked_bin an integer vector keeping track of
#'   which tags belong to a linked bin
#'
#' @slot modal_rearrangement charactering string of the modal
#'   rearrangemnt type
#'
#' @slot percent_rearrangement length-one numeric vector indicating
#'   the percentage of improper read pairs for a linked bin that
#'   belong to the modal rearrangement
#'
#' @slot fraction_linking_tags the fraction of all tags for a linked
#'   bin that link the two tag clusters
#'
#' @slot split_reads a \code{GRanges} object of the split read alignments
#'
#'
#' @export
#' @rdname Rearrangement-class
#' @aliases Rearrangement
setClass("Rearrangement", representation(linkedBins="GRanges",
                                         improper="GAlignmentPairs",
                                         partitioning="integer",
                                         link1id="character",
                                         link2id="character",
                                         tags="GRanges",
                                         tag_map_to_linked_bin="character",
                                         modal_rearrangement="character",
                                         percent_rearrangement="numeric",
                                         fraction_linking_tags="numeric",
                                         split_reads="GRanges"))


conciseGRangeSummary <- function(g, type){
  g <- g[1]
  g2 <- g$linked.to
  chr1 <- chromosome(g)
  chr2 <- chromosome(g2)
  #s1 <- prettyNum(start(g), big.mark=",")
  s2 <- prettyNum(start(g2), big.mark=",")
  e1 <- prettyNum(end(g), big.mark=",")
  #e2 <- prettyNum(end(g2), big.mark=",")
  txt1 <- paste0(chr1, ":", e1, "-", chr2, ":", s2)
  ##txt2 <- paste0("chr2:", s2, "-", e2)
  paste0(substr(type, 1, 1), "||", txt1)
}

setMethod("show", "Rearrangement", function(object){
  cat(paste0("Object of class '", class(object), "'"), "\n")
  cat("  number of rearrangements: ", length(linkedBins(object)), "\n")
  cat("  number of linking RPs: ", length(improper(object)), "\n")
  cat("  number of split reads: ", length(splitReads(object))/2, "\n")
  cat("  number of tags : ", length(tags(object)), "\n")
  if(length(object)==1){
    r <- modalRearrangement(object)
    txt <- conciseGRangeSummary(linkedBins(object), r)
    cat("  ", txt, "\n")
    cat("  modal type   :", r, "\n")
    cat("  %tags for mode :", round(percentRearrangement(object), 0), "\n")
  }
  cat("  see linkedBins(), improper(), and tags() \n")
})




##--------------------------------------------------
##
##
##
## Methods
##
##
##
##--------------------------------------------------

#' @rdname fractionLinkingTags
#' @aliases fractionLinkingTags,Rearrangement-method
setMethod("fractionLinkingTags", "Rearrangement", function(object)
  object@fraction_linking_tags)

#' @rdname fractionLinkingTags
#' @aliases fractionLinkingTags,Rearrangement,ANY-method
setReplaceMethod("fractionLinkingTags", "Rearrangement", function(object, value){
  object@fraction_linking_tags <- value
  object
})

link1id <- function(object) object@link1id
link2id <- function(object) object@link2id

#' @aliases linkedBins,Rearrangement-method
#' @rdname linkedBins-methods
setMethod("linkedBins", "Rearrangement", function(object) object@linkedBins)

#' @aliases, linkedBins,Rearrangement,ANY-method
#' @rdname linkedBins-methods
setReplaceMethod("linkedBins", "Rearrangement", function(object, value){
  object@linkedBins <- value
  ##
  ## list only the improper tags that are within the redefined bin
  ##
  all_tags <- tags(object)
  keep <- overlapsAny(all_tags, value) | overlapsAny(all_tags, value$linked.to)
  tags(object) <- all_tags[keep]
  object
})

#' @rdname modalRearrangement
#' @aliases modalRearrangement,Rearrangement-method
setMethod("modalRearrangement", "Rearrangement", function(object) object@modal_rearrangement)

#' @rdname modalRearrangement
#' @aliases modalRearrangement,Rearrangement,ANY-method
setReplaceMethod("modalRearrangement", "Rearrangement", function(object,value){
  object@modal_rearrangement <- value
  object
})

#' @rdname numberLinkingRP-methods
#' @aliases numberLinkingRP,Rearrangement-method
#' @export
setMethod("numberLinkingRP", "Rearrangement", function(object){
  nsupportingReads <- table(names(improper(object)))
  nsupportingReads <- nsupportingReads[names(object)]
  as.integer(nsupportingReads)
})

#' @aliases partitioning,Rearrangement-method
#' @rdname Rearrangement-class
setMethod("partitioning", "Rearrangement", function(object) object@partitioning)

#' @rdname percentRearrangement
#' @aliases percentRearrangement,Rearrangement-method
setMethod("percentRearrangement", "Rearrangement", function(object) object@percent_rearrangement)

#' @rdname percentRearrangement
#' @aliases percentRearrangement,Rearrangement,ANY-method
setReplaceMethod("percentRearrangement", "Rearrangement", function(object,value){
  object@percent_rearrangement <- value
  object
})


#' @rdname splitReads
#' @aliases splitReads,Rearrangement-method
setMethod("splitReads", "Rearrangement", function(object) object@split_reads)

#' @aliases splitReads,Rearrangement,GRanges-method
#' @rdname splitReads
setReplaceMethod("splitReads", c("Rearrangement", "GRanges"),
                 function(object, value){
                   object@split_reads <- value
                   object
                 })

#' @aliases tags,Rearrangement-method
#' @rdname tags-methods
setMethod("tags", "Rearrangement", function(object) object@tags)

#' @aliases tags,Rearrangement,ANY-method
#' @rdname tags-methods
setReplaceMethod("tags", "Rearrangement", function(object, value){
  object@tags <- value
  lb <- linkedBins(object)
  ## update the tag mapping to the rearrangement bin
  linked.ids <- strsplit(names(lb), "-")
  linked1id <- sapply(linked.ids, "[", 1)
  linked2id <- sapply(linked.ids, "[", 2)
  hits1 <- findOverlaps(value, lb, select="first")
  hits2 <- findOverlaps(value, lb$linked.to, select="first")
  i1 <- linked1id[hits1]
  i2 <- linked2id[hits2]
  mapping <- rep(NA, length(tags))
  mapping[!is.na(i1)] <- i1[!is.na(i1)]
  mapping[is.na(i1)] <- i2[is.na(i1)]
  object@tag_map_to_linked_bin <- mapping
  object
})

setMethod("tagMapping", "Rearrangement", function(object) object@tag_map_to_linked_bin)

#' @rdname Rearrangement-class
#' @aliases length,Rearrangement-method
setMethod("length", "Rearrangement", function(x) length(linkedBins(x)))

#' @rdname Rearrangement-class
#' @aliases names,Rearrangement-method
setMethod("names", "Rearrangement", function(x) names(linkedBins(x)))

##--------------------------------------------------
##
## subsetting
##
##--------------------------------------------------

#' @return a \code{Rearrangement} object
#' @rdname Rearrangement-class
#' @param x a \code{Rearrangement} object
#' @param i a \code{logical} vector
#' @param j ignored
#' @param ... ignored
#' @param drop ignored
setMethod("[", "Rearrangement", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    if(is(i, "array")) i <- as.logical(i)
    lb <- x@linkedBins[i]
    ##partition <- lb$partition
    p <- partitioning(x)
    p <- p[ names(p) %in% names(lb) ]
    irp <- improper(x)
    rpid <- names(irp)
    irp <- irp [ rpid %in% names(lb) ]

    x@linkedBins <- lb
    x@partitioning <- p
    x@improper <- irp

    x@link1id <- x@link1id[i]
    x@link2id <- x@link2id[i]
    index <- tagMapping(x) %in% link1id(x) | tagMapping(x) %in% link2id(x)
    x@tags <- x@tags[ index ]
    x@tag_map_to_linked_bin = tagMapping(x)[ index ]
  }
  x
})
