#' @include ExonSubset-class.R 
NULL

#' Container for storing CDS near a sequence junction
#'
#' @slot tA a \code{GRangesList}
#' @slot tB a \code{GRangesList}
#' @slot tC a \code{GRangesList}
#' @slot tD a \code{GRangesList}
#' @slot id sample id
#' @rdname Transcripts-class
#' @export
setClass("Transcripts", representation(tA="GRangesList",
                                       tB="GRangesList",
                                       tC="GRangesList",
                                       tD="GRangesList",
                                       id="character"))


#' A container for fused transcripts
#'
#' Container for fused transcripts.
#'
#' @examples
#' TranscriptsFusion()
#' @slot fusions a \code{list} of fusions
#' @export
#' @rdname TranscriptsFusion-class
setClass("TranscriptsFusion", representation(fusions="list"),
         contains="Transcripts")

#' A container for clipped transcripts
#'
#' @examples
#' ClippedTranscripts()
#' @slot left a \code{GRangesList}
#' @slot right a \code{GRangesList}
#' @rdname ClippedTranscripts-class
#' @export
setClass("ClippedTranscripts", representation(left="GRangesList", right="GRangesList",
                                              names="character"))

#' @examples
#' Transcripts()
#' @export
#' @keywords-internal
#' @rdname Transcripts-class
Transcripts <- function(tA=GRangesList(),
                        tB=GRangesList(),
                        tC=GRangesList(),
                        tD=GRangesList(),
                        id=character()){
  new("Transcripts", tA=tA, tB=tB, tC=tC, tD=tD, id=id)
}

#' @rdname TranscriptsFusion-class
#' @export
TranscriptsFusion <- function(object, fusions=list()){
  if(missing(object)){
    object <- Transcripts()
  }
  new("TranscriptsFusion",
      tA=txA(object),
      tB=txB(object),
      tC=txC(object),
      tD=txD(object), id=identifier(object), fusions=fusions)
}


#' @aliases ClippedTranscripts,missing,missing-method
#' @rdname ClippedTranscripts-class
#' @export
setMethod("ClippedTranscripts", c("missing", "missing"),
          function(transcripts, i, left=GRangesList(), right=GRangesList(),
                   names=c("5prime", "3prime")){
            new("ClippedTranscripts", left=left, right=right,
                names=names)
          })

setMethod("names", "ClippedTranscripts", function(x) x@names)

setMethod("elementNROWS", "Transcripts", function(x){
  c(length(txA(x)), length(txB(x)), length(txC(x)), length(txD(x)))
})


#' @aliases fusions,TranscriptsFusion-method
#' @rdname TranscriptsFusion-class
#' @export
setMethod("fusions", "TranscriptsFusion", function(object) object@fusions)

#' @aliases fusions,TranscriptsFusion,list-method
#' @rdname TranscriptsFusion-class
#' @export
setReplaceMethod("fusions", c("TranscriptsFusion", "list"), function(object, value){
  object@fusions <- value
  object
})


#' @rdname TranscriptsFusion-class
#' @aliases fusionNames,TranscriptsFusion-method
#' @export
setMethod("fusionNames", "TranscriptsFusion", function(object) names(fusions(object)))

#' @export
#' @rdname Transcripts-class
#' @aliases identifier,Transcripts-method
setMethod("identifier", "Transcripts", function(object) object@id)

#' @export
#' @rdname Transcripts-class
#' @aliases identifier,Transcripts,ANY-methodsetMethod("identifier", "Transcripts", function(object) object@id)


setReplaceMethod("identifier", "Transcripts", function(object, value){
  object@id <- value
  object
})

setMethod("identifier", "Transcripts", function(object) object@id)

setReplaceMethod("identifier", "Transcripts", function(object, value){
  object@id <- value
  object
})

setMethod("identifier", "Transcripts", function(object) object@id)

setReplaceMethod("identifier", "Transcripts", function(object, value){
  object@id <- value
  object
})

setMethod("identifier", "Transcripts", function(object) object@id)

setReplaceMethod("identifier", "Transcripts", function(object, value){
  object@id <- value
  object
})


setReplaceMethod("identifier", "Transcripts", function(object, value){
  object@id <- value
  object
})


#' @aliases left,ClippedTranscripts-method
#' @param x a \code{ClippedTranscripts} object
#' @param ... ignored
#' @rdname ClippedTranscripts-class
setMethod("left", "ClippedTranscripts", function(x, ...) x@left)

#' @aliases numberFusions,TranscriptsFusion-method
#' @rdname TranscriptsFusion-class
#' @export
setMethod("numberFusions", "TranscriptsFusion", function(object) length(fusions(object)))

#' @rdname ClippedTranscripts-class
#' @aliases right,ClippedTranscripts-method
setMethod("right", "ClippedTranscripts", function(x, ...) x@right)

#' @rdname transcript-accessors
#' @aliases txA,Transcripts-method
#' @export
setMethod("txA", "Transcripts", function(object) object@tA)

#' @rdname transcript-accessors
#' @aliases txB,Transcripts-method
#' @export
setMethod("txB", "Transcripts", function(object) object@tB)

#' @rdname transcript-accessors
#' @export
#' @aliases txC,Transcripts-method
setMethod("txC", "Transcripts", function(object) object@tC)

#' @rdname transcript-accessors
#' @export
#' @aliases txD,Transcripts-method
setMethod("txD", "Transcripts", function(object) object@tD)

#' @rdname transcript-accessors
#' @export
#' @aliases txA<-,Transcripts,ANY-method
setReplaceMethod("txA", "Transcripts", function(object, value){
  object@tA <- value
  object
})

#' @rdname transcript-accessors
#' @export
#' @aliases txB<-,Transcripts,ANY-method
setReplaceMethod("txB", "Transcripts", function(object, value){
  object@tB <- value
  object
})

#' @rdname transcript-accessors
#' @export
#' @aliases txC<-,Transcripts,ANY-method
setReplaceMethod("txC", "Transcripts", function(object, value){
  object@tC <- value
  object
})

#' @rdname transcript-accessors
#' @export
#' @aliases txD<-,Transcripts,ANY-method
setReplaceMethod("txD", "Transcripts", function(object, value){
  object@tD <- value
  object
})

#' @rdname Transcripts-class
#' @param x a \code{Transcripts} object
#' @param i integer vector
#' @param j integer vector
#' @param ... ignored
#' @param drop ignored
setMethod("[[", "Transcripts", function(x, i, j, ..., drop=FALSE){
  if(is.character(i)){
    i <- match(i[1], LETTERS[1:4])
  }
  if(i==1) return(txA(x))
  if(i==2) return(txB(x))
  if(i==3) return(txC(x))
  if(i==4) return(txD(x))
})

#' @rdname Transcripts-class
#' @aliases [[<-,Transcripts,ANY,ANY,GRangesList-method
setReplaceMethod("[[", c("Transcripts", "ANY", "ANY", "GRangesList"),
                 function(x, i, j, value){
                   if(is.character(i)){
                     i <- match(i[1], LETTERS[1:4])
                   }
                   if(i==1){
                     x@tA <- value
                   }
                   if(i==2) {
                     x@tB <- value
                   }
                   if(i==3){
                     x@tC <- value
                   }
                   if(i==4){
                     x@tD <- value
                   }
                   return(x)
                 })

setMethod("show", "Transcripts", function(object){
  cat("An object of class 'Transcripts'\n")
  el <- elementNROWS(object)
  cat("   txA: ", el[1], "elements\n")
  cat("   txB: ", el[2], "elements\n")
  cat("   txC: ", el[3], "elements\n")
  cat("   txD: ", el[4], "elements\n")
})

setMethod("show", "TranscriptsFusion", function(object){
  cat("An object of class 'Transcripts'\n")
  el <- elementNROWS(object)
  cat("   txA: ", el[1], "elements\n")
  cat("   txB: ", el[2], "elements\n")
  cat("   txC: ", el[3], "elements\n")
  cat("   txD: ", el[4], "elements\n")
  cat("   number of fusions:", numberFusions(object), "\n")
  cat("   fusion symbols:", paste(names(fusions(object)), collapse=","), "\n")
  cat("   A transcript object can be subset by any of its fusions\n")
  cat("   See fusions()\n")
})

setMethod("show", "ClippedTranscripts", function(object){
  cat("A 'ClippedTranscripts' object\n")
  cat("   clipped on left side of junction:\n", show(left(object)), "\n")
  cat("   clipped on right side of junction:\n", show(right(object)), "\n")
})
