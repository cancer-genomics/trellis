#' @include ExonSubset-class.R 
NULL

#' @aliases ClippedTranscripts
#' @param transcripts a \code{Transcripts} object
#' @param i an integer vector
#' @param left a \code{GRangesList} for the transcripts to the left of the junction
#' @param right a \code{GRangesList} for the transcripts to the right of the junction
#' @rdname ClippedTranscripts-class
#' @export
setGeneric("ClippedTranscripts", function(transcripts, i, left=GRangesList(),
                                          right=GRangesList())
  standardGeneric("ClippedTranscripts"))

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
setClass("ClippedTranscripts", representation(left="GRangesList", right="GRangesList"))

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
          function(transcripts, i, left=GRangesList(), right=GRangesList()){
            new("ClippedTranscripts", left=left, right=right)
          })

anyLengthZero <- function(x)  length(left(x)) == 0 || length(right(x)) == 0

logicalList <- function(x) x@fusions

#' Select CDS involved in a fusion
#'
#' A list of \code{GRangesList} objects, formalized as a
#' \code{Transcripts} class, is subset to return an instance of the
#' same class containing only the CDS believed to be involved in the
#' fusion.
#'
#' @seealso \code{\link{fuse}}
#' 
#' @examples
#'   data(rear_cds)
#'   clip(rear_cds)
#'   ## One of the possible fused transcripts is gene B and gene C, or "BC"
#'   ## There are four possible refseq transcripts associated with gene 'B'
#'   full_refseqs <- txB(rear_cds)
#'   elementNROWS(full_refseqs)
#' 
#'   ## We exclude CDS from the four refseq transcripts that are not
#'   ## involved in the fusion using the clip function:
#'   clipped_tx <- clip(rear_cds)
#'   ## The clipped refseq transcripts for gene B are give by
#'   names(left(clipped_tx))
#'   clipped_refseqs <- left(clipped_tx)
#'   elementNROWS(clipped_refseqs)
#'   ##
#'   ## One transcript on the "right" side
#'   ##
#'   names(right(clipped_tx))
#'   length(right(clipped_tx)[[1]])
#'   ## Note the above exluced 9 CDS that appear in the unclipped refseq transcript
#'   length(txC(rear_cds)[[1]])
#' 
#'   ## The exon ranks of the clipped transcripts shows that the 
#'   ## first 9 CDS are excluded in the fused transcript:
#' 
#'   right(clipped_tx)[[1]]$exon_rank
#' @return A \code{ClippedTranscripts} object
#' @export
#' @param transcripts A \code{TranscriptsFusion} object
#' @param i a length-one integer vector
clip <- function(transcripts, i) {
  logical.list1 <- logicalList(transcripts)[[1]]
  logical.list2 <- logicalList(transcripts)[[2]]
  if(!missing(i)){
    result <- transcripts[logicalList(transcripts)[[i]]]
    return(result)
  }
  nms <- names(fusions(transcripts))
  cl.tx1 <- transcripts[logical.list1]
  cl.tx2 <- transcripts[logical.list2]
  if(anyLengthZero(cl.tx1)){
    return(cl.tx2)
  }
  if(anyLengthZero(cl.tx2)){
    return(cl.tx1)
  }
  stop("Chimeric proteins on both strands. Didn't expect to see this.")
}

setMethod("elementNROWS", "Transcripts", function(x){
  c(length(txA(x)), length(txB(x)), length(txC(x)), length(txD(x)))
})

#' @aliases fuse,TranscriptsFusion-method
#' @rdname fuse-methods
#' @export
setMethod("fuse", "TranscriptsFusion", function(object, nms){
  fuse(clip(object, nms))
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

.getFusedTx <- function(object){
  lt <- left(object)
  rt <- right(object)
  joined <- list()
  it <- 1
  ## if either lt or rt has length zero, joined is an empty list
  for(i in seq_along(lt)){
    L <- lt[[i]]
    nm.lt <- names(lt)[i]
    L$orientation <- "left"
    for(j in seq_along(rt)){
      R <- rt[[j]]
      nm.rt <- names(rt)[j]
      R$orientation <- "right"
      J <- c(L, R)
      joined[[it]] <- J
      names(joined)[it] <- paste(nm.lt, nm.rt, sep="::")
      it <- it+1
    }
  }
  joined <- GRangesList(joined)
}

#' Join two clipped transcripts as a GRangesList object
#'
#' If there are only two transcripts to join, the GRangesList object
#' has length one.  In general, the GRangesList object returned has
#' length R x C where R is the number of transcripts from the gene
#' 5-prime of the fusion and C is the number of transcripts from the
#' gene 3' of the fusion.
#'
#' @seealso \code{link{clip}}
#' @return A \code{GRangesList} 
#' @examples
#'   data(rear_cds)
#'   clipped <- clip(rear_cds)
#'   ## The fused transcript contains exon ranks 2-9 of NM_001288715 and exon
#'   ## ranks 10-57 of NM_007118 
#'   fused <- fuse(clipped)
#' @rdname fuse-methods
#' @export
setMethod("fuse", "ClippedTranscripts", function(object, nms){
  joined <- .getFusedTx(object)
  joined
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


setGeneric("left", function(x, ...) standardGeneric("left"))
setGeneric("right", function(x, ...) standardGeneric("right"))

#' @aliases left,ClippedTranscripts-method
#' @param x a \code{ClippedTranscripts} object
#' @param ... ignored
#' @rdname ClippedTranscripts-class
#' @export
setMethod("left", "ClippedTranscripts", function(x, ...) x@left)

#' @aliases numberFusions,TranscriptsFusion-method
#' @rdname TranscriptsFusion-class
#' @export
setMethod("numberFusions", "TranscriptsFusion", function(object) length(fusions(object)))

#' @rdname ClippedTranscripts-class
#' @aliases right,ClippedTranscripts-method
#'@export
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

## subsetting by one of the fusions produces a ClippedTranscripts object
#' @rdname TranscriptsFusion-class
#' @param x a \code{TranscriptsFusion} object
#' @param i an integer-vector
#' @param j ignored
#' @param ... ignored
#' @param drop ignored
#' @return a \code{ClippedTranscripts} object
setMethod("[", c("TranscriptsFusion", "ExonSubset"), function(x, i, j, ..., drop=FALSE){
  grl1 <- x[[name.left(i)]] ## GRangesList
  ##gr1 <- grl1[[tx1(i)]]  ## GRanges
  grl1 <- grl1[inRearrangement.left(i)]
  grl1 <- grl1[elementNROWS(grl1) > 0]

  grl2 <-x[[name.right(i)]] ## GRangesList
  grl2 <- grl2[inRearrangement.right(i)]  ## GRanges
  grl2 <- grl2[elementNROWS(grl2) > 0]

  x[[name.left(i)]] <- grl1
  x[[name.right(i)]] <- grl2
  ClippedTranscripts(left=grl1, right=grl2)
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
