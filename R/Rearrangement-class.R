## TODO Need examples and tests for constructor

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
#'   which tags belong to which linked bin @slot modal_rearrangement
#'   charactering string of the modal rearrangemnt type
#' 
#' @slot percent_rearrangement length-one numeric vector indicating
#'   the percentage of improper read pairs for a linked bin that
#'   belong to the modal rearrangement
#' 
#' @slot fraction_linking_tags the fraction of all tags for a linked
#'   bin that link the two tag clusters
#' 
#' @export
#' @rdname Rearrangement-class
setClass("Rearrangement", representation(linkedBins="GRanges",
                                         improper="GAlignmentPairs",
                                         partitioning="integer",
                                         link1id="character",
                                         link2id="character",
                                         tags="GRanges",
                                         tag_map_to_linked_bin="character",
                                         modal_rearrangement="character",
                                         percent_rearrangement="numeric",
                                         fraction_linking_tags="numeric"))



##--------------------------------------------------
##
## generics
##
##--------------------------------------------------

#' Accessor for improper read pairs
#'
#' Extract the improperly paired reads.  Reads are flagged as improper
#' by the alignment algorithm.  Typically, this indicates that the
#' separation between a read and its mate is larger than expected, the
#' orientation of the read and its mate is different from that
#' anticipated in the reference genome, or the reads align to
#' different chromosomes.
#' 
#' 
#' @export
#' @rdname improper-methods
#' @seealso \code{\linkS4class{Rearrangement}}
#' @param object a \code{Rearrangement} or \code{StructuralVariant} object
setGeneric("improper", function(object) standardGeneric("improper"))

#' Accessor for clusters of reads belonging to improper read pairs that are linked.
#'
#' Clusters of improper reads that are linked by the read pairs are
#' suggestive of a new (somatic) sequence junction not present in the
#' reference genome.  The linking is represented by a \code{GRanges}
#' object with an interval demarcating the boundaries of a cluster of
#' improper tags.  The tag cluster that is linked by the pairing
#' information is put in the element metadata (\code{mcols}) of the
#' \code{GRanges} object.
#' 
#' @export
#' @rdname linkedBins-methods
#' @param object a \code{Rearrangement} object
setGeneric("linkedBins", function(object) standardGeneric("linkedBins"))

#' @export
#' @rdname linkedBins-methods
#' @aliases linkedBins<-
#' @param value a \code{GRanges} object with a column \code{linked.to}
#'   in the element metadata.
setGeneric("linkedBins<-", function(object,value) standardGeneric("linkedBins<-"))

#' @export
#' @rdname Rearrangement-class
#' @keywords internal
setGeneric("partitioning", function(object) standardGeneric("partitioning"))

#' Accessor for all the improper reads
#'
#' Accessor for improper reads belonging to one of the linked
#' bins. While a given read can only belong to one cluster (the
#' clusters are non-overlapping), a single cluster can be linked to
#' several other clusters.  Consequently, a single read (ignoring its
#' mate) can belong to several linked clusters. For visualizing the
#' reads supporting a rearrangement, it may also be important to carry
#' forward tags that do not support the rearrangement.
#' 
#' @export
#' 
#' @param object a \code{Rearrangement} object
#' @rdname tags-methods
setGeneric("tags", function(object) standardGeneric("tags"))

setGeneric("tagMapping", function(object) standardGeneric("tagMapping"))


link1id <- function(object) object@link1id
link2id <- function(object) object@link2id

#' @aliases improper,Rearrangement-method
#' @rdname improper-methods
setMethod("improper", "Rearrangement", function(object) object@improper)

#' @aliases linkedBins,Rearrangement-method
#' @rdname linkedBins-methods
setMethod("linkedBins", "Rearrangement", function(object) object@linkedBins)

#' @aliases partitioning,Rearrangement-method
#' @rdname Rearrangement-class
setMethod("partitioning", "Rearrangement", function(object) object@partitioning)

#' @aliases tags,Rearrangement-method
#' @rdname tags-methods
setMethod("tags", "Rearrangement", function(object) object@tags)

setMethod("tagMapping", "Rearrangement", function(object) object@tag_map_to_linked_bin)

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


         
