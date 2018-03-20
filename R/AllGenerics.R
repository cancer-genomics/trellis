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


