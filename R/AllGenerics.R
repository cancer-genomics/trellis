#' @include help.R
NULL

##
## Transcripts-class generics
##

#' @aliases ClippedTranscripts
#' @param transcripts a \code{Transcripts} object
#' @param i an integer vector
#' @param left a \code{GRangesList} for the transcripts to the left of the junction
#' @param right a \code{GRangesList} for the transcripts to the right of the junction
#' @rdname ClippedTranscripts-class
#' @export
setGeneric("ClippedTranscripts", function(transcripts, i, left=GRangesList(),
                                          right=GRangesList(), names=c("5prime", "3prime"))
  standardGeneric("ClippedTranscripts"))

#' Fuse transcripts
#'
#' @export
#' @rdname fuse-methods
#' @aliases fuse
#' @param object a \code{Transcripts}-derived class
#' @param nms name of transcript (character-string)
setGeneric("fuse", function(object, nms) standardGeneric("fuse"))

#' @rdname TranscriptsFusion-class
#' @export
setGeneric("fusions", function(object) standardGeneric("fusions"))

#' @rdname TranscriptsFusion-class
#' @param value 
#' @export
setGeneric("fusions<-", function(object, value) standardGeneric("fusions<-"))

#' @rdname TranscriptsFusion-class
#' @export
setGeneric("fusionNames", function(object) standardGeneric("fusionNames"))

#' @rdname Transcripts-class
#' @keywords internal
#' @export
setGeneric("identifier", function(object) standardGeneric("identifier"))

#' @rdname TranscriptsFusion-class
#' @keywords internal
#' @param value character string
#' @export
setGeneric("identifier<-", function(object, value) standardGeneric("identifier<-"))

#' @rdname TranscriptsFusion-class
#' @export
setGeneric("numberFusions", function(object) standardGeneric("numberFusions"))

#' Accessor for transcript from Transcripts-derived class
#'
#' @examples
#' data(rear_cds)
#' txA(rear_cds)
#' txB(rear_cds)
#' txC(rear_cds)
#' txD(rear_cds)
#' @aliases txA txB txC txD
#' @rdname transcript-accessors
#' @export
#' @param object a \code{Transcripts}-derived object
setGeneric("txA", function(object) standardGeneric("txA"))

#' @rdname transcript-accessors
#' @export
setGeneric("txB", function(object) standardGeneric("txB"))

#' @rdname transcript-accessors
#' @export
setGeneric("txC", function(object) standardGeneric("txC"))

#' @rdname transcript-accessors
#' @export
setGeneric("txD", function(object) standardGeneric("txD"))

#' @rdname transcript-accessors
#' @param value a \code{GRangesList} 
#' @export
setGeneric("txA<-", function(object,value) standardGeneric("txA<-"))

#' @rdname transcript-accessors
#' @export
setGeneric("txB<-", function(object,value) standardGeneric("txB<-"))

#' @rdname transcript-accessors
#' @export
setGeneric("txC<-", function(object,value) standardGeneric("txC<-"))

#' @rdname transcript-accessors
#' @export
setGeneric("txD<-", function(object,value) standardGeneric("txD<-"))

setGeneric("left", function(x, ...) standardGeneric("left"))

setGeneric("right", function(x, ...) standardGeneric("right"))

##
## Amplicon-class generics
##
#' @export
#' @rdname AmpliconGraph-class
setGeneric("amplicons", function(object) standardGeneric("amplicons"))

#' @export
#' @rdname AmpliconGraph-class
setGeneric("amplicons<-", function(object, value) standardGeneric("amplicons<-"))

#' @export
#' @rdname AmpliconGraph-class
setGeneric("ampliconRanges", function(object) standardGeneric("ampliconRanges"))

#' @export
#' @rdname AmpliconGraph-class
setGeneric("ampliconRanges<-", function(object,value) standardGeneric("ampliconRanges<-"))

#' @rdname AmpliconGraph-class
#' @export
#' @param object a \code{AmpliconGraph} object
setGeneric("assemblyGaps", function(object) standardGeneric("assemblyGaps"))

#' @rdname AmpliconGraph-class
#' @export
setGeneric("driver", function(object) standardGeneric("driver"))

#' @rdname AmpliconGraph-class
#' @export
setGeneric("graph", function(object) standardGeneric("graph"))

#' @rdname AmpliconGraph-class
#' @export
setGeneric("graph<-", function(object, value) standardGeneric("graph<-"))

#' @rdname AmpliconGraph-class
#' @export
setGeneric("isAmplicon", function(object) standardGeneric("isAmplicon"))

#' @rdname AmpliconGraph-class
#' @export
setGeneric("isSomatic", function(object) standardGeneric("isSomatic"))

#' @rdname AmpliconGraph-class
#' @export
setGeneric("nonampliconRanges", function(object) standardGeneric("nonampliconRanges"))

#' @rdname AmpliconGraph-class
#' @export
setGeneric("queryRanges", function(object) standardGeneric("queryRanges"))

#' @rdname AmpliconGraph-class
#' @export
setGeneric("queryRanges<-", function(object, value) standardGeneric("queryRanges<-"))

##
## Deletion-class Generics
##

#' Type of deletion
#'
#' An accessor for the deletion calls.  The deletion calls are
#' hemizygous, hemizygous+, homozygous, and homozygous+.  The '+'
#' denotes that improper read pairs support the deletion identified
#' from the read depth analysis.
#'
#' @export
#' @rdname calls-method
#' @param object a \code{StructuralVariant}
setGeneric("calls", function(object) standardGeneric("calls"))

#' @rdname calls-method
#' @export
#' @param value a character string indicating the type of deletion
setGeneric("calls<-", function(object, value) standardGeneric("calls<-"))


#' Accessor for variable indicating a grouping for deletions
#'
#' The grouping can indicate, for example, overlapping hemizygous
#' deletions.
#' 
#' @rdname groupedVariant-method
#' @export
#' @param object a \code{StructuralVariant}
#' @return factor indicating a grouping for deletions
#' @keywords internal
setGeneric("groupedVariant", function(object) standardGeneric("groupedVariant"))

#' @rdname groupedVariant-method
#' @export
#' @param value a factor indicating a grouping for deletions
#' @keywords internal
setGeneric("groupedVariant<-", function(object,value) standardGeneric("groupedVariant<-"))

#' Accessor of an index for the proper read pairs
#' 
#' @export
#' @rdname indexing-methods
#' @keywords internal
#' @param object a \code{StructuralVariant}
setGeneric("indexProper", function(object) standardGeneric("indexProper"))

#' @rdname indexing-methods
#' @param value a list
#' @export
setGeneric("indexProper<-", function(object, value) standardGeneric("indexProper<-"))

#' @keywords internal
#' @rdname indexing-methods
#' @export
setGeneric("indexImproper<-", function(object,value) standardGeneric("indexImproper<-"))

#' @return a numeric vector
#' @keywords internal
#' @rdname indexing-methods
#' @export
setGeneric("indexImproper", function(object) standardGeneric("indexImproper"))

#' Get / set properly paired reads
#'
#' @return \code{GAlignmentPairs}
#' @export
#' @keywords internal
#' @rdname proper-methods
setGeneric("proper", function(object) standardGeneric("proper"))

#' @rdname proper-methods
#' @param value a \code{GAlignmentPairs} object
#' @export
setGeneric("proper<-", function(object, value) standardGeneric("proper<-"))

#' @rdname StructuralVariant-class
#' @keywords internal
#' @export
setGeneric("readPairs", function(object) standardGeneric("readPairs"))

#' Extract deletion regions from a StructuralVariant object
#'
#' @return a \code{GRanges} object for the genomic intervals of the deletions
#' 
#' @export
#' @param object a \code{StructuralVariant}
#' @rdname variant-methods
setGeneric("variant", function(object) standardGeneric("variant"))

#' @param value a \code{GRanges} object
#' @rdname variant-methods
#' @export
setGeneric("variant<-", function(object,value) standardGeneric("variant<-"))


##--------------------------------------------------
##
## Rearrangement-class generics
##
##--------------------------------------------------

#' The fraction of improper read pairs that link two clusters of improper reads
#'
#' For bona fide deletions, inversions, and translocations, we would
#' expect that the fraction of improper read pairs that link two
#' clusters of improper reads to be high.  For amplicons, the fraction
#' of reads linking two clusters may be small.
#' 
#' @return numeric
#' @export
#' @param object a \code{Rearrangement} object
setGeneric("fractionLinkingTags", function(object) standardGeneric("fractionLinkingTags"))

#' @rdname fractionLinkingTags
#' @param value length-one numeric vector
#' @export
setGeneric("fractionLinkingTags<-", function(object,value)
  standardGeneric("fractionLinkingTags<-"))

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

#' Accessor for the modal rearrangement of a linked tag cluster
#'
#' For a two clusters of improper reads that are linked by the pairing
#' information, we classify the type of rearrangement supported by
#' each improper read pair.  Once each read pair is typed, we store
#' the modal type in the \code{Rearrangement} object.
#'
#' 
#' @return a character string
#' @export
#' @param object a \code{Rearrangement} object
setGeneric("modalRearrangement", function(object) standardGeneric("modalRearrangement"))
#' @export

#' @rdname modalRearrangement
#' @param value a character-string indicating the modal rearrangement type
#' @export
setGeneric("modalRearrangement<-", function(object, value) standardGeneric("modalRearrangement<-"))

#' Calculate the number of read pairs that link two clusters of improper reads
#'
#' Clusters of improper read pairs can be linked by the mates of the
#' improper reads that belong to these clusters.  This function
#' computes the number of improperly paired reads that link two
#' clusters. The number of read pairs that would link two clusters in
#' a bona fide somatic rearrangement depends on the depth of
#' sequencing.  For 30x coverage, we require at least 5 linking read
#' pairs.
#' 
#' @rdname numberLinkingRP-methods
#' @keywords internal
#' @return a numeric vector
#' @export
setGeneric("numberLinkingRP", function(object) standardGeneric("numberLinkingRP"))

#' @export
#' @rdname Rearrangement-class
#' @keywords internal
setGeneric("partitioning", function(object) standardGeneric("partitioning"))

#' The fraction of improper read pairs supporting the modal rearrangement type
#'
#' @seealso \code{\link{modalRearrangement}}
#' @export
#' @param object a \code{Rearrangement} object
setGeneric("percentRearrangement", function(object) standardGeneric("percentRearrangement"))

#' @rdname percentRearrangement
#' 
#' @param value a numeric vector with range [0,1] indicating the
#'   fraction of improper read pairs that support the modal
#'   rearrangement
#' @export
setGeneric("percentRearrangement<-",
           function(object, value) standardGeneric("percentRearrangement<-"))


#' Accessor for split read alignments
#'
#' @rdname splitReads
#'
#' @param object a \code{Rearrangement} object
#'
#' @export
setGeneric("splitReads",
           function(object) standardGeneric("splitReads"))

#' @rdname splitReads
#'
#' @param value a \code{GRanges} object of split read alignments
#'
#' @export
setGeneric("splitReads<-",
           function(object, value) standardGeneric("splitReads<-"))


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

#' @rdname tags-methods
#' @param value a \code{GRanges} object of read alignments
#' @keywords internal
#' @export
setGeneric("tags<-", function(object, value) standardGeneric("tags<-"))

setGeneric("tagMapping", function(object) standardGeneric("tagMapping"))

setGeneric("first<-", function(x, value) standardGeneric("first<-"))

setGeneric("last<-", function(x, value) standardGeneric("last<-"))

##
## DeletionList-class Generics
##

#' GRangesList-derived class for deletions
#'
#' Each element in the list are the set of deletions for a sample
#' 
#' @export
#' @rdname DeletionList-methods
#' @param object see \code{showMethods("DeletionList")}
setGeneric("DeletionList", function(object) standardGeneric("DeletionList"))


##
## ExonSubset-class Generics
##

#' @export
#' @keywords-internal
#' @rdname ExonSubset-class
#' @param object a \code{ExonSubset} object
#' @return a \code{LogicalList}
setGeneric("inRearrangement.left", function(object) standardGeneric("inRearrangement.left"))

#' @export
#' @keywords-internal
#' @rdname ExonSubset-class
#' @return a \code{LogicalList}
setGeneric("inRearrangement.right", function(object) standardGeneric("inRearrangement.right"))

#' @export
#' @keywords-internal
#' @rdname ExonSubset-class
setGeneric("name.left", function(object) standardGeneric("name.left"))

#' @export
#' @keywords-internal
#' @rdname ExonSubset-class
setGeneric("name.right", function(object) standardGeneric("name.right"))
