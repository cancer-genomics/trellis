#' @include help.R
NULL

#' Parameter class for rearrangement analysis
#'
#' @slot rp_separation length-one numeric vector
#' @slot min_number_tags_per_cluster length-one numeric vector
#' @slot min_cluster_size length-one numeric vector
#' @slot max_cluster_size length-one numeric vector
#' @slot min.gapwidth length-one numeric vector
#' @slot percent_modal_type length-one numeric vector
#' @slot percent_linking length-one numeric vector
#' @export
#' @rdname RearrangementParams-class
setClass("RearrangementParams", representation(rp_separation="numeric",
                                               min_number_tags_per_cluster="numeric",
                                               min_cluster_size="numeric",
                                               max_cluster_size="numeric",
                                               min.gapwidth="numeric",
                                               percent_modal_type="numeric",
                                               percent_linking="numeric"))

#' A parameter class for somatic rearrangement analysis
#'
#' Some details about why we look for linked tag clusters.
#'
#' @details A tag cluster is defined as follows:
#'
#' (i)   it must have at least <min_number_tags_per_cluster> reads
#' 
#' (ii) Provided (i) is TRUE, the cluster includes all improper reads
#'   with less than <min.gapwidth> separation
#'
#' (iii) The size of a cluster is defined as the difference in the
#'   minimum basepair across all members and the maximum basepair
#'   across all members.  The size of the cluster must be at least
#'   <min_cluster_size> and no bigger than <max_cluster_size>.
#'
#' Having determined the type of rearrangment supported by each read
#'   pair for two linked clusters, we require that the modal
#'   rearrangement type be supported by at least <percent_modal_type>>
#'   read pairs.
#'
#'
#' @examples
#' 
#' ## Default rearrangement parameters for whole genome sequence data
#' ##  with 30x coverage
#' rp <- RearrangementParams()
#' 
#' @export
#' 
#' @param rp_separation length-one numeric vector indicating minimum
#'   separation of the first and last read of a pair
#' 
#' @param min_number_tags_per_cluster length-one numeric vector
#'   indicating the minimum number of reads in a cluster
#' 
#' @param min_cluster_size length-one numeric vector; the minimum size
#'   of a cluster of reads
#' 
#' @param max_cluster_size length-one numeric vector; the maximum size
#'   of a cluster of reads
#' 
#' @param min.gapwidth length-one numeric vector; reads with at most
#'   min.gapwidth separation between them are considered overlapping.
#' 
#' @param percent_modal_type length-one numeric vector; the percentage
#'   of reads that must agree with the modal rearrangement type
#'
#' @param percent_linking length-one numeric vector; two linked tag
#'   clusters must be linked by at least this percentage of reads. See
#'   details
#'
#' @rdname RearrangementParams-class
RearrangementParams <- function(rp_separation=10e3,
                                min_number_tags_per_cluster=5,
                                min_cluster_size=115L,
                                max_cluster_size=5000L,
                                min.gapwidth=1000L,
                                percent_modal_type=0.9,
                                percent_linking=0.8){
  new("RearrangementParams", rp_separation=rp_separation,
      min_number_tags_per_cluster=min_number_tags_per_cluster,
      min_cluster_size=min_cluster_size,
      max_cluster_size=max_cluster_size,
      min.gapwidth=min.gapwidth,
      percent_modal_type=percent_modal_type,
      percent_linking=percent_linking)
}

#' @param object a \code{RearrangementParams} object
#' @rdname RearrangementParams-class
#' @export
#' @keywords internal
setGeneric("rpSeparation", function(object) standardGeneric("rpSeparation"))

#' @rdname RearrangementParams-class
#' @export
#' @keywords internal
setGeneric("minNumberTagsPerCluster", function(object) standardGeneric("minNumberTagsPerCluster"))

#' @rdname RearrangementParams-class
#' @export
#' @keywords internal
setGeneric("minClusterSize", function(object) standardGeneric("minClusterSize"))

#' @rdname RearrangementParams-class
#' @export
#' @keywords internal
setGeneric("maxClusterSize", function(object) standardGeneric("maxClusterSize"))

#' @rdname RearrangementParams-class
#' @export
#' @keywords internal
setGeneric("minGapWidth", function(object) standardGeneric("minGapWidth"))

#' @rdname RearrangementParams-class
#' @export
setGeneric("percentModalType", function(object) standardGeneric("percentModalType"))

#' @rdname RearrangementParams-class
#' @aliases rpSeparation,RearrangementParams-method
setMethod("rpSeparation", "RearrangementParams", function(object) object@rp_separation)

#' @rdname RearrangementParams-class
#' @aliases minNumberTagsPerCluster,RearrangementParams-method
setMethod("minNumberTagsPerCluster", "RearrangementParams", function(object) object@min_number_tags_per_cluster)

#' @rdname RearrangementParams-class
#' @aliases minNumberTagsPerCluster,RearrangementParams-method
setMethod("minClusterSize", "RearrangementParams", function(object) object@min_cluster_size)

#' @rdname RearrangementParams-class
#' @aliases maxClusterSize,RearrangementParams-method
setMethod("maxClusterSize", "RearrangementParams", function(object) object@max_cluster_size)

#' @rdname RearrangementParams-class
#' @aliases minGapWidth,RearrangementParams-method
setMethod("minGapWidth", "RearrangementParams", function(object) object@min.gapwidth)

#' @rdname RearrangementParams-class
#' @aliases percentModalType,RearrangementParams-method
setMethod("percentModalType", "RearrangementParams", function(object) object@percent_modal_type)

#' @rdname RearrangementParams-class
#' @export
percentLinking <- function(object) object@percent_linking

setMethod("show", "RearrangementParams", function(object){
  cat("Object of class RearrangementParams \n")
  cat("  min tag separation:", rpSeparation(object), "\n")
  cat("  min tags/cluster:", minNumberTagsPerCluster(object), "\n")
  cat("  min cluster size:", minClusterSize(object), "\n")
  cat("  max cluster size:", maxClusterSize(object), "\n")
  cat("  min gap width between tags:", minGapWidth(object), "\n")
  cat("  prop modal rearrangement:", percentModalType(object), "\n")
  cat("  prop linking 2 clusters :", percentLinking(object), "\n")
})
