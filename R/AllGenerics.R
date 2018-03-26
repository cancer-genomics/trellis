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
#' @param names character vector indicating strand of transcripts ('5prime' or '3prime').
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
#' @param value a \code{GRanges} object
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

#' @export
#' @param x a \code{Transcripts}-derived object
#' @param ... ignored
#' @rdname transcript-accessors
setGeneric("left", function(x, ...) standardGeneric("left"))

#' @export
#' @rdname transcript-accessors
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

#' @rdname improper-methods
#' @export
setGeneric("improper<-", function(object, value) standardGeneric("improper<-"))

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

## import from S4Vectors
##setGeneric("first<-", function(x, value) standardGeneric("first<-"))

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


#' Set numeric scale for \code{PreprocessViews2} object
#'
#' Exported only for internal use by other packages.
#' 
#' @return numeric
#' 
#' @export
#' @docType methods
#' @param x A \code{PreprocessViews2} object
#' @param value a length-one numeric vector
#' @rdname setScale-method
#' @keywords internal
setGeneric("setScale<-", function(x, value) standardGeneric("setScale<-"))

##
## Preprocess-class generics
##

#' Indexes the GRanges of a BamViews-derived class
#'
#' @docType methods
#' @rdname indexRanges-method
#' @param object a \code{BamViews}-derived object
#' @export
#' @keywords internal
setGeneric("indexRanges", function(object) standardGeneric("indexRanges"))

#' @param value an integer-vector 
#' @rdname indexRanges-method
#' @export
#' @keywords internal
setGeneric("indexRanges<-", function(object, value) standardGeneric("indexRanges<-"))

#' Accessor for file paths
#'
#' File path to intermediate data saved to disk
#' 
#' @export
#' @param object a \code{PreprocessViews2} object
#' @rdname paths
#' @aliases paths paths<-
setGeneric("paths", function(object) standardGeneric("paths"))

#' @param value a character-vector of file paths to intermediate data
#' @rdname paths
#' @aliases paths,PreprocessViews2-method  paths<-,PreprocessViews2-method
#' @export
setGeneric("paths<-", function(object,value) standardGeneric("paths<-"))


##
## RearrangementParams-class generics
##

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

##
## RearrangementList-class generics
##

#' Constructor for \code{RearrangementList} class
#'
#' @return a \code{RearrangementList} object
#' @rdname RearrangementList-class
#' @export
#' @keywords internal
#' @param object a list of \code{Rearrangement} objects
#'
#' @param modal_rearrangement a character vector of the modal
#'   rearrangement corresponding to each element in \code{object}.
#'   Must be the same length as \code{object}.
#'
#' @param percent_rearrangement a numeric vector indicating the
#'   fraction of improper read pairs supporting the modal
#'   rearrangement type.  Must be the same length as \code{object}.
#'
#' @param colData a \code{DataFrame} containing metadata on the
#'   \code{Rearrangement} objects.
#'
setGeneric("RearrangementList", function(object,
                                         modal_rearrangement,
                                         percent_rearrangement,
                                         colData) standardGeneric("RearrangementList"))


##
## SegmentParam-class generics
##

#' @param x \code{SegmentParam} object
#' @rdname SegmentParam-class
#' @export
setGeneric("cbs_alpha", function(x) standardGeneric("cbs_alpha"))

#' @rdname SegmentParam-class
#' @export 
setGeneric("undo.splits", function(x) standardGeneric("undo.splits"))

#' @rdname SegmentParam-class
#' @export 
setGeneric("undo.SD", function(x) standardGeneric("undo.SD"))


##
## copy number generics
##

#' Accessor for 'copy' assays
#'
#' Extract matrix of log2-transformed estimates of copy number
#' relative to autosomal mode or median
#'
#' @export
#' @docType methods
#' @param object A \code{RangedSummarizedExperiment}
#' @rdname copynumber-methods
setGeneric("copynumber", function(object) standardGeneric("copynumber"))

#' @rdname copynumber-methods
#' @export
setGeneric("copynumber<-", function(object, value) standardGeneric("copynumber<-"))

##
## miscellaneous generics
##

#' Extract autosomal (human) sequence levels
#' 
#' @export
#' @param object a \code{GRanges} object
#' @docType methods
#' @rdname chromosome-methods
setGeneric("autosomeNames", function(object) standardGeneric("autosomeNames"))

#' Return seqnames of a \code{GRanges} object as a character vector
#' 
#' @export
#' @docType methods
#' @rdname chromosome-methods
setGeneric("chromosome", function(object) standardGeneric("chromosome"))

setGeneric("modev", function(x, na.rm=TRUE) standardGeneric("modev"))

##
## Generics for CNVs  (svcnvs)
##
setGeneric("zeroEndAnchors", function(object) standardGeneric("zeroEndAnchors"))
setGeneric("singleEndAnchors", function(object) standardGeneric("singleEndAnchors"))
setGeneric("bothEndAnchors", function(object) standardGeneric("bothEndAnchors"))

#' Remove the genomic intervals in query that overlap with intervals
#' in subject
#'
#'
#' @param query a \code{\linkS4class{BamViews}} or \code{\linkS4class{GRanges}} instance
#' @param subject a \code{\linkS4class{GRanges}} object
#' @param type see \code{\link{findOverlaps}}
#' @param ... Additional arguments passed to \code{findOverlaps}
#' @return Returns the \code{query} with intervals that overlap \code{subject} removed
#' @export
#' @seealso \code{\link[GenomicRanges]{findOverlaps}}
#' @rdname filterBy-methods
setGeneric("filterBy", function(query, subject, type="any", ...) standardGeneric("filterBy"))

setGeneric("clusterReadPairs", function(object) standardGeneric("clusterReadPairs"))
setGeneric("reviseJunction", function(object) standardGeneric("reviseJunction"))
setGeneric("variant<-", function(object, value) standardGeneric("variant<-"))
setGeneric("CNAObject", function(object, valuename) standardGeneric("CNAObject"))

##
## Generics from svalignments
##

#' Extract sequences of reads supporting rearrangements from a bam file
#'
#'
#' @export
#' @return a data.frame
#' @param object a \code{Rearrangement} or \code{RearrangementList} object
#' @param bam.file complete path to BAM file
#' @param params a \code{RearrangementParams} object
#'
#' @param MAX the maximum number of read pairs to extract for a
#'   rearrangement.  If the number of read pairs supporting a
#'   rearrangement is greater than MAX, a random sample of MAX
#'   supporting read pairs is returned.
#'
#' @param build the reference genome buld that reads were aligned to.  Currently
#' supported builds include "hg19" and "hg18".
#'
#' @rdname getSequenceOfReads-methods
setGeneric("getSequenceOfReads", function(object, bam.file,
                                          params=RearrangementParams(),
                                          MAX=25L, build)
  standardGeneric("getSequenceOfReads"))


# #' Parse BAM file for improper read pairs near a set of GRanges
# #'
# #' All reads aligned to the intervals given by
# #' \code{queryRanges(object)} are identified by the low-level function
# #' \code{.scan_all_readpairs}.  This function reads alignments by
# #' \code{readGAlignments} and then makes pairs of the alignments by
# #' \code{makeGAlignmentPairs2}.  The latter function is an adaption of
# #' the function \code{makeGAlignmentPairs} implemented in the
# #' \code{GenomeAlignments} package but allows for the read pairs to be
# #' improper.
# #'
# #' @param object Typically an \code{AmpliconGraph}, though the only
# #'   requirement is that the method \code{queryRanges} is defined
# #' @param bam.file character-vector providing valid complete path to a
# #'   bam file
# #' @param flags length-two integer vector as given by \code{scanBamFlags}
# #' @export
# setGeneric(name = "get_readpairs", def = function(object, bam.file, flags=scanBamFlag()) {standardGeneric("get_readpairs")})
# setMethod(
#   f = "get_readpairs",
#   signature = "AmpliconGraph",
#   definition = function(object, bam.file, flags=scanBamFlag()) {
#     g <- queryRanges(object)
#     .get_readpairs2(g, bam.file, flags)
#   }
# )
# setMethod(
#   f = "get_readpairs",
#   signature = "GRanges",
#   definition = function(object, bam.file, flags=scanBamFlag()) {
#     .get_readpairs2(g, bam.file, flags)
#   }
                                        # )

setGeneric("rename", function(x, ...) standardGeneric("rename"))

#' Accessor for the combination of strands involved in a rearrangement
#'
#' @rdname Rearrangement-class
#' @aliases type
#' @export
setGeneric("type", function(object) standardGeneric("type"))
setGeneric("last<-", function(x, value) standardGeneric("last<-"))
setGeneric("numberLinkingRP", function(object) standardGeneric("numberLinkingRP"))

#' A region linked by improperly paired reads
#'
#' @param x a \code{GRanges} object with value `linked.to` in the `mcols`
#' @export
setGeneric("linkedTo", function(x) standardGeneric("linkedTo"))
