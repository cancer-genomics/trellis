#' @include AllGenerics.R
NULL

getExonsAmp2_4 <- function(left, right, transcripts){
  fusion.list <- list(ExonSubset(name.left="C",
                                 inRearrangement.left=end(txC(transcripts)) <
                                   start(right),
                                 name.right="A",
                                 inRearrangement.right=start(txA(transcripts)) >
                                   end(left)),
                      ExonSubset(name.left="B",
                                 inRearrangement.left=start(txB(transcripts)) >
                                   end(left),
                                 name.right="D",
                                 inRearrangement.right=end(txD(transcripts)) <
                                   start(right)))
  names(fusion.list) <- sapply(fusion.list, names)
  tx <- TranscriptsFusion(transcripts, fusions=fusion.list)
  tx
}

getExonsDel1_2 <- function(left, right, transcripts){
  fusion.list <- list(ExonSubset(name.left="A",
                                 inRearrangement.left=start(txA(transcripts)) <
                                   start(left),
                                 name.right="C",
                                 inRearrangement.right=start(txC(transcripts)) >
                                   end(right)),
                      ExonSubset(name.left="D",
                                 inRearrangement.left=start(txD(transcripts)) >
                                   end(right),
                                 name.right="B",
                                 inRearrangement.right=start(txB(transcripts)) <
                                   start(left)))
  names(fusion.list) <- sapply(fusion.list, names)
  tx <- TranscriptsFusion(transcripts, fusions=fusion.list)
  tx
}

getExonsInversionBalanced <- function(left, right, transcripts){
  fusion.list <- list(ExonSubset(name.left="A",
                                 inRearrangement.left=end(txA(transcripts)) <
                                   start(left),
                                 name.right="D",
                                 inRearrangement.right=end(txD(transcripts)) <
                                   start(right)),
                      ExonSubset(name.left="C",
                                 inRearrangement.left=start(txC(transcripts)) <
                                   start(right),
                                 name.right="B",
                                 inRearrangement.right=end(txB(transcripts)) <
                                   start(left)),
                      ExonSubset(name.left="D",
                                 inRearrangement.left=end(txD(transcripts)) >
                                   end(right),
                                 name.right="A",
                                 inRearrangement.right=start(txA(transcripts)) >
                                   end(left)),
                      ExonSubset(name.left="B",
                                 inRearrangement.left=start(txB(transcripts)) >
                                   end(left),
                                 name.right="C",
                                 inRearrangement.right=end(txC(transcripts)) >
                                   end(right)))
  names(fusion.list) <- sapply(fusion.list, names)
  tx <- TranscriptsFusion(transcripts, fusions=fusion.list)
  tx
}

getExonsInversion1_2 <- function(left, right, transcripts){
  ##
  ## C and D flip. Translating 4 pieces :  presumes a stop after the
  ## last D exon and the last A exon in the rearranged genome.
  ##
  fusion.list <- list(ExonSubset(name.left="A",
                                 inRearrangement.left=end(txA(transcripts)) < start(left),
                                 name.right="D",
                                 inRearrangement.right=end(txD(transcripts)) < start(right)),
                      ExonSubset(name.left="C",
                                 inRearrangement.left=end(txC(transcripts)) < start(right),
                                 name.right="B",
                                 inRearrangement.right=end(txB(transcripts)) < start(left)))
  names(fusion.list) <- sapply(fusion.list, names)
  tx <- TranscriptsFusion(transcripts, fusions=fusion.list)
  tx
}

getExonsInversion3_4 <- function(left, right, transcripts){
  fusion.list <- list(ExonSubset(name.left="B",
                                 inRearrangement.left=start(txB(transcripts)) > end(left),
                                 name.right="C",
                                 inRearrangement.right=start(txC(transcripts)) > end(right)),
                      ExonSubset(name.left="D",
                                 inRearrangement.left=start(txD(transcripts)) > end(right),
                                 name.right="A",
                                 inRearrangement.right=start(txA(transcripts)) > end(left)))
  names(fusion.list) <- sapply(fusion.list, names)
  tx <- TranscriptsFusion(transcripts, fusions=fusion.list)
  tx
}

getExonsTranslocationBalanced <- function(left, right, transcripts){
  ## if there is a balanced translocation and the tag clusters for the
  ## breakpoints overlap, all 4 types are compatible
  ## if this is the case, we should see R1+ reads on 2 different
  ## chromosomes
  fusion.list <- list(ExonSubset(name.left="A",
                                 inRearrangement.left=end(txA(transcripts)) < start(left),
                                 name.right="C",
                                 inRearrangement.right=start(txC(transcripts)) > end(right)),
                      ExonSubset(name.left="D",
                                 inRearrangement.left=start(txD(transcripts)) > end(right),
                                 name.right="B",
                                 inRearrangement.right=end(txB(transcripts)) < start(left)),
                      ExonSubset(name.left="C",
                                 inRearrangement.left=end(txC(transcripts)) < start(right),
                                 name.right="A",
                                 inRearrangement.right=start(txA(transcripts)) > end(left)),
                      ExonSubset(name.left="B",
                                 inRearrangement.left=end(txB(transcripts)) > end(left),
                                 name.right="D",
                                 inRearrangement.right=end(txD(transcripts)) < start(right)))
  names(fusion.list) <- sapply(fusion.list, names)
  tx <- TranscriptsFusion(transcripts, fusions=fusion.list)
  tx
}


getCDFun <- function(call){
  parsecall <- gsub(",", "_", call)
  fun <- switch(parsecall,
                deletion1_deletion2=getExonsDel1_2,
                amp2_amp4=getExonsAmp2_4,
                trans1_trans2_trans3_trans4=getExonsTranslocationBalanced,
                inv1_inv2=getExonsInversion1_2,
                inv3_inv4=getExonsInversion3_4,
                inv5_inv6_inv7_inv8=getExonsInversionBalanced,
                NULL)
  if(is.null(fun)) stop("call not recognized")
  fun
}

charstrand <- function(x) as.character(strand(x))

grlapply <- function(X, FUN, ...) GRangesList(lapply(X, FUN, ...))

rearrangementTranscripts <- function(gr, tx, cds){
  left <- gr
  right <- gr$linked.to
  tx1 <- subsetByOverlaps(tx, left, ignore.strand=TRUE, maxgap=5000)
  tx2 <- subsetByOverlaps(tx, right, ignore.strand=TRUE, maxgap=5000)
  i1 <- match(as.character(tx1$tx_name), names(cds))
  i1 <- i1[!is.na(i1)]
  i2 <- match(as.character(tx2$tx_name), names(cds))
  i2 <- i2[!is.na(i2)]
  if(length(i1) > 0){
    cds.left <- cds[i1]
  } else cds.left <- GRanges()
  if(length(i2) > 0){
    cds.right <- cds[i2]
  } else cds.right <- GRanges()
  cdsA <- grlapply(cds.left, function(x){
    x <- x[charstrand(x) == "+"]
  })
  if(all(elementNROWS(cdsA)==0))
    cdsA <- GRangesList()
  cdsB <- grlapply(cds.left, function(x) {
    x <- x[charstrand(x) == "-"]
  })
  if(all(elementNROWS(cdsB)==0))
    cdsB <- GRangesList()
  cdsC <- grlapply(cds.right, function(x){
    x <- x[charstrand(x) == "+"]
  })
  if(all(elementNROWS(cdsC)==0))
    cdsC <- GRangesList()
  cdsD <- grlapply(cds.right, function(x){
    x <- x[charstrand(x) == "-"]
  })
  if(all(elementNROWS(cdsD)==0))
    cdsD <- GRangesList()
  ta <- cdsA[elementNROWS(cdsA) > 0]
  tb <- cdsB[elementNROWS(cdsB) > 0]
  tc <- cdsC[elementNROWS(cdsC) > 0]
  td <- cdsD[elementNROWS(cdsD) > 0]
  Transcripts(tA=ta, tB=tb, tC=tc, tD=td)
}

rearrangementTranscripts2 <- function(gr, tx, cds){
  left <- gr
  right <- gr$linked.to
  tx1 <- subsetByOverlaps(tx, left, ignore.strand=TRUE, maxgap=5000)
  tx2 <- subsetByOverlaps(tx, right, ignore.strand=TRUE, maxgap=5000)
  i1 <- match(as.character(tx1$tx_name), names(cds))
  i1 <- i1[!is.na(i1)]
  i2 <- match(as.character(tx2$tx_name), names(cds))
  i2 <- i2[!is.na(i2)]
  if(length(i1) > 0){
    cds.left <- cds[i1]
  } else cds.left <- GRanges()
  if(length(i2) > 0){
    cds.right <- cds[i2]
  } else cds.right <- GRanges()
  list("5p"=cds.left, "3p"=cds.right)
##  cdsA <- grlapply(cds.left, function(x){
##    x <- x[charstrand(x) == "+"]
##  })
##  if(all(elementNROWS(cdsA)==0))
##    cdsA <- GRangesList()
##  cdsB <- grlapply(cds.left, function(x) {
##    x <- x[charstrand(x) == "-"]
##  })
##  if(all(elementNROWS(cdsB)==0))
##    cdsB <- GRangesList()
##  cdsC <- grlapply(cds.right, function(x){
##    x <- x[charstrand(x) == "+"]
##  })
##  if(all(elementNROWS(cdsC)==0))
##    cdsC <- GRangesList()
##  cdsD <- grlapply(cds.right, function(x){
##    x <- x[charstrand(x) == "-"]
##  })
##  if(all(elementNROWS(cdsD)==0))
##    cdsD <- GRangesList()
##  ta <- cdsA[elementNROWS(cdsA) > 0]
##  tb <- cdsB[elementNROWS(cdsB) > 0]
##  tc <- cdsC[elementNROWS(cdsC) > 0]
##  td <- cdsD[elementNROWS(cdsD) > 0]
##  Transcripts(tA=ta, tB=tb, tC=tc, tD=td)
}

#' Finds all genes overlapping a rearrangement (ignoring strand) and
#' returns the full CDS of each
#'
#' Finds all genes overlapping the \code{linkedBins} of a
#' \code{Rearrangement} object (ignoring strand) and extracts the full
#' CDS.  The rearrangement type inferred from the position and strand
#' orientation of read pairs (see \code{RearrangementType}) determines
#' the orientation of the genes in the rearranged genome.  See detail
#' for the generic symbols used to denote a fused transcript.
#'
#' @details
#'
#' #' \preformatted{
#' Reference   (-- chr1, == chr2,  .... not in rearranged transcript)
#'           A                                    C
#' 5' + ------------|.........   ...............|==============
#' 3' - ------------|.........   ...............|==============
#'           B                                    D
#'
#'
#' Tumor
#'
#'       chr1  A        C
#' 5' + -------------|==========
#' 3' - -------------|==========
#'             B        D
#'
#' }
#'
#' \strong{Deletions:} 
#' For a deletion, we have chr1 = chr2, R1+ < R2-, and R1- > R2+.  The
#'   possible fusions are AC and DB.
#'
#' \strong{Translocations:} A read pair in a translocation supports
#'   one of \code{trans1}, \code{trans2}, \code{trans3}, and
#'   \code{trans4}.  For a balanced translocation, all four types
#'   could be observed.  We analyze the data as a balanced
#'   translocation and evaluate all possible fusions even if not all
#'   types are observed.  Again, the possible fusions are AC and DB.
#'
#' Amplicons and translocations in which the sequences is inverted in
#'   the rearranged genome are analyzed as inversions.
#'
#' \strong{Inversions:} There are 8 possible fusions indicated by an
#'   inverted read pair denoted as \code{invX}, where X = 1-8.
#'   Possible fused gene products are AD, CB, BC, and DA.
#'
#'
#'
#' @examples
#'   options(warn=-1)
#'   library("TxDb.Hsapiens.UCSC.hg19.refGene")
#'   txdb <- TxDb.Hsapiens.UCSC.hg19.refGene
#'   tx <- transcripts(txdb)
#'   cds.all <- cdsBy(txdb, "tx", use.names=TRUE)
#'   data(rear_list, package="trellis")
#'   r <- rear_list[["18557-18736"]]
#'   ##data(rear_cds, package="trellis")
#'   ##trace(getCDS, browser)
#'   rear_cds <- getCDS(r, tx, cds.all)
#'   fusions(rear_cds)
#'   ## There are for transcripts for gene 'B'.  The accessor txB
#'   ## returns the CDSs for each of the transcripts as a
#'   ## \code{GRangesList}
#'   txB(rear_cds)
#'   ## Similarly, there are accessors \code{txA}, \code{txC}, and
#'   ## \code{txD}
#' @return a \code{Transcripts} object
#' @export
#' @param rear A \code{Rearrangement} object
#' @param tx A \code{GRanges} object of transcripts from the \code{TxDb} object
#' @param cds A \code{GRangesList} object containing the CDS for each transcript
getCDS <- function(rear, tx, cds){
  tx <- keepSeqlevels(tx, seqlevels(improper(rear)), pruning.mode="coarse")
  trx <- rearrangementTranscripts(linkedBins(rear), tx, cds)
  identifier(trx) <- names(linkedBins(rear))
  if(all(elementNROWS(trx)==0)) return(trx)
  ##
  ## here 'left' and 'right' refer to the left tag cluster and the
  ## right tag cluster.  If a large gene spans both tag clusters, then
  ## left and right would both include the complete list of exons for
  ## this gene.  Typically, a gene will only span one tag cluster and
  ## either left will be empty or right will be empty.
  ##
  ## In particular, left and right does not imply the inferred
  ## orientation after the rearrangement
  left.cluster <- granges(linkedBins(rear))
  right.cluster <- granges(linkedBins(rear)$linked.to)
  modal.call <- modalRearrangement(rear)
  cdfun <- getCDFun(modal.call)
  cds.object <- cdfun(left.cluster, right.cluster, trx)
  cds.object
}

sequentialExonRanks <- function(cds){
  cds <- endoapply(cds, function(x) {
    x$original_exon_rank <- x$exon_rank
    x$exon_rank <- seq_along(x)
    x
  })
  cds
}

.trim_one_base <- function(dna.ss){
  i <- length(dna.ss)
  dna.ss <- dna.ss[-i]
  dna.ss
}

dnaSeqLengths <- function(dna.ss) sapply(dna.ss, length)

codonMultiple <- function(dna.ss){
  dnaSeqLengths(dna.ss) %% 3 == 0
}

endSeqAtCodon <- function(dna.ss){
  index <- which(!codonMultiple(dna.ss))
  dna.ss[index] <- endoapply(dna.ss[index], .trim_one_base)
  dna.ss
}

divisibleBy3Seq <- function(dna.ss){
  tmp <- dna.ss
  if(all(codonMultiple(dna.ss))) return(dna.ss)
  dna.ss <- endSeqAtCodon(dna.ss)
  if(all(codonMultiple(dna.ss))) return(dna.ss)  
  dna.ss <- endSeqAtCodon(dna.ss)
  ## should only have to apply endSeqAtCodon twice
  if(!all(codonMultiple(dna.ss))) {
    stop()
  }
  dna.ss
}

extractTranscriptSeqs2 <- function(genome, fused.txlist){
  seqs <- getSeq(genome, fused.txlist)
  seqs <- DNAStringSet(lapply(seqs, unlist))
  seqs <- divisibleBy3Seq(seqs)
  seqs
}

#' Translate the sequence of a rearranged DNA sequence
#'
#' @examples
#' data(rear_cds, package="trellis")
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' genome <- BSgenome.Hsapiens.UCSC.hg19
#' fused_tx <- fuse(clip(rear_cds))
#' fused.proteins <- tumorProtein(genome, fused_tx)
#'
#' @param genome a \code{BSgenome} object
#' @param fused.txlist a \code{GRangesList}
#' @export
#' @seealso \code{\link{referenceProtein}}
tumorProtein <- function(genome, fused.txlist){
  ## Extract the DNA sequence of the fused transcripts
  fused.txlist2 <- sequentialExonRanks(fused.txlist)
  seqlevels(fused.txlist2, pruning.mode="coarse") <- seqlevelsInUse(fused.txlist2)
  if(length(seqlevels(fused.txlist2)) > 1){
    fused.seq <- extractTranscriptSeqs2(genome, fused.txlist2)
  } else {
    fused.seq <- extractTranscriptSeqs2(genome, fused.txlist2)
  }
  fused.protein <- translate(fused.seq, if.fuzzy.codon="solve")
}

#' Translate the unrearranged (reference) sequence
#'
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' genome <- BSgenome.Hsapiens.UCSC.hg19
#' data(rear_cds, package="trellis")
#' unrear_cds <- fullTranscripts(rear_cds)
#' ##
#' ## referenceProtein
#' ##
#' ref.proteins <- referenceProtein(genome, unrear_cds, names(unrear_cds))
#' @param genome a \code{BSgenome} object
#' @param cds a \code{GRangesList} of transcript CDS
#' @param nms names of transcripts to translate
#' @export
#' @seealso \code{\link{tumorProtein}}
referenceProtein <- function(genome, cds, nms){
  ## Extract the DNA sequence of the original transcripts
  ## Reference transcripts
  nms <- gsub("\\(promoter\\)", "", nms)
##  any.promoter <- any(nms == "promoter")
##  if(any.promoter){
##    nms <- nms[nms != "promoter"]
##  }
  ix <- match(nms, names(cds))
  ref.txlist <- cds[ix]
  ref.seq <- extractTranscriptSeqs2(genome, ref.txlist)
  ref.protein <- translate(ref.seq, if.fuzzy.codon="solve")
  ref.protein
}

#' List full transcripts for genes spanning a rearrangement
#'
#' @param object a \code{Transcripts} object
#' @return a \code{GRangesList}
#' @examples
#' data(rear_cds, package="trellis")
#' fullTranscripts(rear_cds)
#' @export
#' @seealso \code{\link{inFrameFusions}}
fullTranscripts <- function(object){
  cds <- list(txA(object),
              txB(object),
              txC(object),
              txD(object))
  cds2 <- do.call("c", cds)
  cds2
}

numberLeftCodons <- function(tx.list){
  ncodons.left <- sapply(tx.list, function(x){
    x <- x[x$orientation == "left"]
    floor(sum(width(x))/3)
  })
  ncodons.left
}

rightCodonIndices <- function(tx.list){
  ncodons.left <- numberLeftCodons(tx.list)
  ncodons.total <-  floor(sum(width(tx.list))/3)
  codons.right <- vector("list", length(ncodons.left))
  for(i in seq_along(ncodons.left)){
    codons.right[[i]] <- seq(ncodons.left[i]+1, ncodons.total[i])
  }
  names(codons.right) <- names(tx.list)
  codons.right
}

isInFrame <- function(query, ref){
  query <- as.character(query)
  ref <- as.character(ref)
  grepl(query, ref, fixed=TRUE)
}



.isInFrame <- function(tumor.protein, ref.protein, nm, i){
  id <- strsplit(nm, "::")[[1]][[2]]
  isInFrame(tumor.protein[i], ref.protein[[id]])
}

#' Evaluates whether a rearranged DNA sequence is in-frame
#'
#' A rearranged sequence is considered in-frame if the protein sequence derived
#' from the DNA sequence 3-prime of the new sequence junction is a subset of the
#' reference protein sequence.
#' 
#' @examples
#'   library(BSgenome.Hsapiens.UCSC.hg19)
#'   genome <- BSgenome.Hsapiens.UCSC.hg19
#'   data(rear_cds, package="trellis")
#'   unrear_cds <- fullTranscripts(rear_cds)
#'   ##
#'   ## referenceProtein
#'   ##
#'   ref_protein <- referenceProtein(genome, unrear_cds, names(unrear_cds))
#'   fused_tx <- fuse(clip(rear_cds))
#'   fused_protein <- tumorProtein(genome, fused_tx)
#'   inFrameFusions(fused_protein, ref_protein, fused_tx)
#' @export
#' @param fused.proteins the rearranged protein sequence
#' @param ref.protein the unrearranged (reference) protein sequence
#' @param fused.txlist a \code{GRangesList} of the fused transcripts
#'
#' @seealso See \code{\link{referenceProtein}} and \code{\link{tumorProtein}}
#'   for extracting the unrearranged and rearranged protein sequences,
#'   respectively. See \code{\link{clip}} and
#'   \code{\link{fuse}} for extracting the fused, clipped transcript
#'   sequence as a \code{GRangesList} object. See \code{\link{fullTranscripts}}
#'   for extracting the full transcripts as a \code{GRangesList}.
inFrameFusions <- function(fused.proteins, ref.protein, fused.txlist){
  codons.right <- rightCodonIndices(fused.txlist)
  is_in_frame <- vector("list", length(fused.proteins))
  for(j in seq_along(fused.proteins)){
    fused.protein <- fused.proteins[[j]]
    nm <- names(fused.proteins)[j]
    i <- codons.right[[j]]
    is_in_frame[[j]] <- .isInFrame(fused.protein, ref.protein, nm, i)
  }
  unlist(is_in_frame)
}

anyPromoter <- function(x){
  xx <- unique(unlist(strsplit(x, "::")))
  length(grep("promoter", xx)) > 0
}

promoterTxName <- function(x){
  xx <- unique(unlist(strsplit(x, "::")))
  xx <- xx[[1]]
  tx <- strsplit(xx, "_")[[1]][2]
  tx
}

.fusionTable <- function(fused.txlist, fused.proteins,
                         in_frame, org.db, txdb,
                         linkedbin.id, id=""){
  tx.nms1 <- unique(sapply(strsplit(names(fused.proteins), "::"), "[", 1))
  tx.nms1 <- unique(gsub("\\(promoter\\)", "", tx.nms1))
  tx.nms2 <- unique(sapply(strsplit(names(fused.proteins), "::"), "[", 2))
  tx.nms2 <- unique(gsub("\\(promoter\\)", "", tx.nms2))
  tx.nms <- unique(c(tx.nms1, tx.nms2))
  i <- grep("promoter", names(fused.txlist))
  is.promoter <- rep(FALSE, length(fused.txlist))
  is.promoter[i] <- TRUE
  ##tx.nms <- unique(gsub("\\(promoter\\)", "", tx.nms))
  ##  any.promoter <- any(tx.nms == "promoter")
  ##  if(any.promoter){
  ##    tx.nms <- tx.nms[tx.nms != "promoter"]
  ##  }
  refseq.map <- suppressMessages(select(org.db, keys=tx.nms,
                                        columns=c("GENENAME",
                                                  "SYMBOL", "REFSEQ",
                                                  "MAP"),
                                        keytype="REFSEQ"))
  rownames(refseq.map) <- tx.nms
  tx1 <- tx.nms1
  tx2 <- tx.nms2
  gene1 <- refseq.map[tx1, ]
  gene2 <- refseq.map[tx2, ]
  if(any(is.na(gene2$SYMBOL))) browser()
  map <- refseq.map$MAP
  chr <- paste0("chr", sapply(strsplit(map, "[pq]"), "[", 1))
  chr.gene1 <- chr[1]
  chr.gene2 <- chr[2]
  ##chr.gene1 <- refseq.map$TXCHROM[1]
  ##chr.gene2 <- refseq.map$TXCHROM[2]
  exons.gene1 <- rep(NA, length(fused.txlist))
  cdsStart.gene1 <- exons.gene1
  cdsEnd.gene1 <- exons.gene1
  strand.gene1 <- exons.gene1
  exons.gene1[!is.promoter] <- unlist(sapply(fused.txlist[!is.promoter], function(x){
    x <- x[x$orientation=="left"]
    paste(min(x$exon_rank), max(x$exon_rank), sep="-")
  }))
  cdsStart.gene1[!is.promoter] <- unlist(sapply(fused.txlist[!is.promoter], function(x){
    x <- x[x$orientation=="left"]
    min(start(x))
  }))
  cdsEnd.gene1[!is.promoter] <- unlist(sapply(fused.txlist[!is.promoter], function(x){
    x <- x[x$orientation=="left"]
    max(end(x))
  }))
  strand.gene1[!is.promoter] <- unlist(sapply(fused.txlist[!is.promoter], function(x){
    x <- x[x$orientation=="left"]
    unique(as.character(strand(x)))
  }))
  exons.gene2 <- sapply(fused.txlist, function(x){
    x <- x[x$orientation=="right"]
    paste(min(x$exon_rank), max(x$exon_rank), sep="-")
  })
  cdsStart.gene2 <- sapply(fused.txlist, function(x){
    x <- x[x$orientation=="right"]
    min(start(x))
  })
  cdsEnd.gene2 <- sapply(fused.txlist, function(x){
    x <- x[x$orientation=="right"]
    max(end(x))
  })
  strand.gene2 <- sapply(fused.txlist, function(x){
    x <- x[x$orientation=="right"]
    unique(as.character(strand(x)))
  })
  ##if(any(is.na(gene2$SYMBOL))) browser()
  tab <- data.frame(fusion=names(fused.proteins),
                    inframe=in_frame,
                    gene1=gene1$SYMBOL,
                    gene1name=gene1$GENENAME,
                    gene2=gene2$SYMBOL,
                    gene2name=gene2$GENENAME,
                    chr.gene1=chr.gene1,
                    chr.gene2=chr.gene2,
                    exonrank.gene1=exons.gene1,
                    exonrank.gene2=exons.gene2,
                    cdsStart.gene1=cdsStart.gene1,
                    cdsEnd.gene1=cdsEnd.gene1,
                    cdsStart.gene2=cdsStart.gene2,
                    cdsEnd.gene2=cdsEnd.gene2,
                    rearrangement.id=linkedbin.id,
                    strand.gene1=strand.gene1,
                    strand.gene2=strand.gene2,
                    stringsAsFactors=FALSE,
                    id=id)
  rownames(tab) <- NULL
  tab
}

#' Determine all possible fusions of a Rearrangement object
#'
#' Determine all possible fusions of a \code{Rearrangement} object and
#' evaluate whether each is in-frame.
#'
#' @details Fusions are analyzed as follows.  First, all CDS from
#'   genes overlapping a rearrangement are extracted (ignoring strand)
#'   using the \code{getCDS} function.  In addition to extracting all
#'   the CDS, this function returns the possible gene fusion that may
#'   result from a rearrangement based on the modal rearrangement type
#'   (the modal rearrangement type is inferred from the strand and
#'   position orientation of read pairs).  The genes are denoted
#'   generically by
#'
#' \preformatted{
#'             A          C
#'  5+ --------------|---------
#'  3- --------B-----|----D----
#'
#' where "|" denotes a new sequence junction in a rearranged genome.
#'   } For a given fusion (say AC), we then \code{clip} CDS from A and
#'   CDS from B that are absent in the fused product.  After clipping,
#'   we \code{fuse} the remaining CDS from genes A and C. The function
#'   \code{tumorProtein} is used to derive the amino acid sequence of
#'   the tumor protein -- the protein that would be formed by fusing
#'   the tripped CDS from genes A and C.  To assess whether the fusion
#'   is in frame, we extract all known full transcripts from genes A
#'   and C and translate the DNA sequence of each transcript to an
#'   amino acid sequence. We refer to the amino acid sequences of the
#'   full CDS as the reference protein.  The function
#'   \code{referenceProtein} is a wrapper for getting the reference
#'   amino acid sequences.  Given the amino acid sequence of the
#'   clipped and fused transcripts (fused tumor protein) and the amino
#'   acid sequence of the full, unclipped transcripts (reference
#'   protein), we compare their sequences to assess whether the fusion
#'   is in-frame using the function \code{inFrameFusions}. The results
#'   are summarized in tabular format by the function
#'   \code{.fusionTable}.
#' 
#' 
#'
#' @seealso See \code{\link{getCDS}} for how the CDS from genes
#'   involved in a rearrangement are extracted, \code{\link{clip}} and
#'   \code{\link{fuse}} for how transcripts are clipped and then
#'   fused, respectively.  See \code{\link{referenceProtein}} and
#'   \code{\link{tumorProtein}} for deriving germline (unrearranged)
#'   and somatic (rearranged) amino acid sequences.
#' 
#' @examples
#' library(org.Hs.eg.db)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' library(TxDb.Hsapiens.UCSC.hg19.refGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.refGene
#' genome <- BSgenome.Hsapiens.UCSC.hg19
#' tx <- transcripts(txdb)
#' options(warn=-1)
#' cds.all <- cdsBy(txdb, "tx", use.names=TRUE)
#' data(rear_list)
#' r <- rear_list[["18557-18736"]]
#' fusionTable(r, txdb, tx, cds.all, genome,
#'             org.Hs.eg.db, id="test")
#' ## in-frame fusion of CTNND2 and TRIO
#' @export
#' @param robj A \code{Rearrangement} object
#' @param txdb A \code{TxDb} object
#' @param tx A \code{GRanges} object of transcripts from the \code{TxDb} object
#' @param cds A \code{GRangesList} object containing the CDS for each transcript
#' @param genome A \code{BSgenome} object
#' @param orgdb A \code{OrdDb} object for the Hsapiens genome
#' @param id A length-one character vector of the sample identifier
fusionTable <- function(robj, txdb, tx, cds, genome, orgdb, id=""){
  cds.fusions <- getCDS(robj, tx, cds)
  is.fusion <- isFusion(cds.fusions)
  if(!is.fusion){
    return(NULL)
  }
  if(sum(elementNROWS(cds.fusions)) <= 1) return(NULL)
  if(is.null(clip(cds.fusions))) return(NULL)
  fused.txlist <- fuse(clip(cds.fusions))
  if(length(fused.txlist) == 0) return(NULL)
  fused.proteins <- tumorProtein(genome, fused.txlist)
  tx.nms <- unique(unlist(strsplit(names(fused.proteins), "::")))
  cds <- fullTranscripts(cds.fusions)
  ref.proteins <- referenceProtein(genome, cds, tx.nms)
  in_frame <- inFrameFusions(fused.proteins, ref.proteins, fused.txlist)
  tab <- .fusionTable(fused.txlist=fused.txlist,
                      fused.proteins=fused.proteins,
                      in_frame=in_frame,
                      org.db=orgdb,
                      txdb=txdb,
                      linkedbin.id=names(robj),
                      id=id)
  tab
}

loadOrgDb <- function(){
  orgdb.path <- system.file("extdata", package="org.Hs.eg.db")
  orgdb.file <- file.path(orgdb.path, "org.Hs.eg.sqlite")
  orgdb <- loadDb(orgdb.file, packageName="org.Hs.eg.db")
  orgdb
}

overlapsAnyTranscript <- function(rlist, build=c("hg19", "hg18"), ...){
  build <- match.arg(build)
  txdb <- loadTxDb(build)
  tx <- transcripts(txdb)
  ## NR_ identifiers are noncoding RNA
  tx <- tx[grep("NM_", tx$tx_name)]
  overlapsAny(linkedBins(rlist), tx, ...) |
    overlapsAny(linkedTo(rlist), tx, ...)
}


## .ftable <- function(robj, id=""){
##   cds.fusions <- getCDS(robj, tx, cds)
##   is.fusion <- isFusion(cds.fusions)
##   if(!is.fusion){
##     return(NULL)
##   }
##   if(sum(elementNROWS(cds.fusions)) <= 1) return(NULL)
##   if(is.null(clip(cds.fusions))) return(NULL)
##   fused.txlist <- fuse(clip(cds.fusions))
##   if(length(fused.txlist) == 0) return(NULL)
##   fused.proteins <- tumorProtein(genome, fused.txlist)
##   tx.nms <- unique(unlist(strsplit(names(fused.proteins), "::")))
##   cds <- fullTranscripts(cds.fusions)
##   ref.proteins <- referenceProtein(genome, cds, tx.nms)
##   in_frame <- inFrameFusions(fused.proteins, ref.proteins, fused.txlist)
##   tab <- .fusionTable(fused.txlist=fused.txlist,
##                       fused.proteins=fused.proteins,
##                       in_frame=in_frame,
##                       org.db=orgdb,
##                       txdb=txdb,
##                       linkedbin.id=names(robj),
##                       id=id)
##   tab
## }

## try using with lexical scope
fTable2 <- function(robj, build, id){
  txdb <- loadTxDb(build)
  orgdb <- loadOrgDb()
  bs.pkg <- paste0("BSgenome.Hsapiens.UCSC.", build)
  genome <- getBSgenome(bs.pkg)
  tx <- transcripts(txdb)
  cds <- suppressWarnings(cdsBy(txdb, "tx", use.names=TRUE))
  getCDS <- getCDS
  rearrangementTranscripts <- rearrangementTranscripts
  tumorProtein <- tumorProtein
  referenceProtein <- referenceProtein
  ftable <- fusionTable
  ftable(robj, id)
}


#' A title
#'
#' @param rlist a \code{RearrangementList}
#' @param txdb a \code{TxDb} object
#' @param orgdb a \code{OrgDb} object
#' @export
fusionList <- function(rlist, id="", build=c("hg19", "hg18")){
  build <- match.arg(build)
  txdb <- loadTxDb(build)
  orgdb <- loadOrgDb()
  bs.pkg <- paste0("BSgenome.Hsapiens.UCSC.", build)
  genome <- getBSgenome(bs.pkg)
  tx <- transcripts(txdb)
  rlist <- rlist[overlapsAnyTranscript(rlist, build, maxgap=5000)]
  if(length(rlist) > 0){
    tx <- keepSeqlevels(tx, seqlevels(rlist), pruning.mode="coarse")
  } else{
    tx <- keepSeqlevels(tx, paste0("chr", c(1:22, "X", "Y")),
                        pruning.mode="coarse")
  }
  cds <- suppressWarnings(cdsBy(txdb, "tx", use.names=TRUE))
  tab.list <- vector("list", length(rlist))
  for(j in seq_along(rlist)){
    tab <- fusionTable(rlist[[j]],
                       txdb=txdb,
                       tx=tx,
                       cds=cds,
                       genome=genome,
                       orgdb=orgdb,
                       id=id)
    tab.list[[j]] <- tab
  }
  tab.list <- tab.list[!sapply(tab.list, is.null)]
  tab <- do.call("rbind", tab.list)
  tab
}

## #' Evaluates rearrangements for in-frame fusions
## #'
## #' Saves the results of fusionTable to disk.  If file already exists,
## #' only read in the file and return a list of \code{tbl_df} objects.
## #'
## #' @seealso \code{fusionTable}
## #'
## #' @export
## #' @return a list of \code{tbl_df} objects
## #' @param file file path
## #' @param rlist a \code{RearrangementList}
## #' @param dirs a \code{DataPaths} object
## #' @param txdb a \code{TxDb} object
## #' @param orgdb an \code{OrgDb} object
## #' @param genome aa \code{BSgenome} object
## #' @param id sample id (character string)
## fusionExperiment <- function(file, id, rlist, dirs, txdb, orgdb, genome) {
##   if(length(rlist)==0){
##     return(NULL)
##   }
##   if(missing(txdb)) stop("missing txdb")
##   if(file.exists(file)){
##     z <- readRDS(file)
##     if(!is.null(z)){
##       z <- tbl_df(z)
##     }
##     return(z)
##   }
##   tab <- fusionList(rlist, txdb=txdb,
##                     genome=genome,
##                     orgdb=orgdb,
##                     id=id)
##   saveRDS(tab, file=file)
##   tbl_df(tab)
## }

anyLengthZero <- function(x)  length(left(x)) == 0 || length(right(x)) == 0

bothLengthZero <- function(x)  length(left(x)) == 0 && length(right(x)) == 0

logicalList <- function(x) x@fusions

#' Select CDS involved in a fusion
#'
#' A list of \code{GRangesList} objects, formalized as a \code{Transcripts}
#' class, is subset to return an instance of the same class containing only the
#' CDS believed to be involved in the fusion.
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
clip <- function(transcripts) {
  logical.list1 <- logicalList(transcripts)[[1]]
  logical.list2 <- logicalList(transcripts)[[2]]
  lens1 <- lengths(logical.list1)
  lens2 <- lengths(logical.list2)
  if(all(lens1) > 0){
    cl.tx <- transcripts[logical.list1]
  }
  if(all(lens2) > 0){
    cl.tx <- transcripts[logical.list2]
  }
  if(!exists("cl.tx")) stop("no clipped transcript generated")
  return(cl.tx)
##  nms <- names(fusions(transcripts))
##  cl.tx1 <- transcripts[logical.list1]
##  cl.tx2 <- transcripts[logical.list2]
##
##  if(bothLengthZero(cl.tx1)){
##    return(cl.tx2)
##  }
##  if(bothLengthZero(cl.tx2)){
##    return(cl.tx1)
##  }
##  stop("Chimeric proteins on both strands. Didn't expect to see this.")
}

#' @aliases fuse,TranscriptsFusion-method
#' @rdname fuse-methods
#' @export
setMethod("fuse", "TranscriptsFusion", function(object, nms){
  fuse(clip(object), nms)
})

.getFusedTx <- function(object){
  lt <- left(object)
  rt <- right(object)
  joined <- list()
  it <- 1
  ## if either lt or rt has length zero, joined is an empty list
  ## RS: 12/13/2016  Why should joined be an empty list??
  ## e.g., the promoter of the left transcript fused with the right transcript
  ##if(length(lt) == 0){
  if(all(elementNROWS(lt)==0)){
    for(j in seq_along(rt)){
      ##
      ## the promoter could be for any of the transcripts -- we'll just use the
      ## name of the first
      ##
      nm.lt <- paste0(names(lt[1]), "(promoter)")
      R <- rt[[j]]
      nm.rt <- names(rt)[j]
      R$orientation <- "right"
      ##J <- c(L, R)
      J <- R
      joined[[it]] <- J
      names(joined)[it] <- paste(nm.lt, nm.rt, sep="::")
      it <- it+1
    }
    joined <- GRangesList(joined)
    return(joined)
  }
  for(i in seq_along(lt)){
    L <- lt[[i]]
    nm.lt <- names(lt)[i]
    if(length(L) > 0){
      L$orientation <- "left"
    } else{
      nm.lt <- paste0(nm.lt, "(promoter)")
    }
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

isFusion <- function(tx.fusion){
  x <- fusions(tx.fusion)
  x1 <- x[[1]]
  lt1 <- inRearrangement.left(x1)
  rt1 <- inRearrangement.right(x1)

  x2 <- x[[2]]
  lt2 <- inRearrangement.left(x2)
  rt2 <- inRearrangement.right(x2)

  (length(lt1) > 0 && length(rt1) > 0) ||
    (length(lt2) > 0 && length(rt2) > 0)
}

#' Subsetting by one of the fusions produces a ClippedTranscripts object
#'
#' @rdname TranscriptsFusion-class
#' @param x a \code{TranscriptsFusion} object
#' @param i an integer-vector
#' @param j ignored
#' @param ... ignored
#' @param drop ignored
#' @return a \code{ClippedTranscripts} object
setMethod("[", c("TranscriptsFusion", "ExonSubset"), function(x, i, j, ..., drop=FALSE){
  name.left <- name.left
  name.right <- name.right
  grl1 <- x[[name.left(i)]] ## GRangesList
  ##gr1 <- grl1[[tx1(i)]]  ## GRanges
  grl1 <- grl1[inRearrangement.left(i)]
  ## this way we keep the names of promoters
  ##grl1 <- grl1[elementNROWS(grl1) > 0]

  grl2 <-x[[name.right(i)]] ## GRangesList
  grl2 <- grl2[inRearrangement.right(i)]  ## GRanges
  grl2 <- grl2[elementNROWS(grl2) > 0]
  x[[name.left(i)]] <- grl1
  x[[name.right(i)]] <- grl2
  ClippedTranscripts(left=grl1, right=grl2)
})

setMethod("lengths", "ExonSubset", function(x, use.names=TRUE) {
  c(length(inRearrangement.left(x)), length(inRearrangement.right(x)))
})

cstrand <- function(x) as.character(strand(x))
isNegStrand <- function(x) cstrand(x)[1] == "-"

use_negative_coords <- function(g){
  g2 <- GRanges(chromosome(g), IRanges(-1*end(g), -1*start(g)))
  g2
}


clip5prime <- function(tx, jxn){
  if(isNegStrand(tx)){
    tx2 <- use_negative_coords(tx)
    jxn["5p"] <- use_negative_coords(jxn["5p"])
  } else tx2 <- tx
  select <- end(tx2) < start(jxn["5p"])
  tx[select]
}


clipFivePrime <- function(tx, jxn){
  if(isNegStrand(tx)){
    tx2 <- use_negative_coords(tx)
    jxn <- use_negative_coords(jxn)
    select <- end(tx2) < start(jxn)
    return(tx[select])
  }
  select <- end(tx) < start(jxn)
  tx[select]
}

clipFivePrimeJunctions <- function(cds, coding_jxns){
  tx5p.nms <- fivePrimeTranscriptNames(coding_jxns)
  rid.all <- names(tx5p.nms)
  cds5p.all <- cds[ tx5p.nms ]
  cds5p.list <- split(cds5p.all, rid.all)
  cds5p.truncated <- vector("list", length(coding_jxns))
  names(cds5p.truncated) <- names(cds5p.list)
  for(i in seq_along(cds5p.list)){
    rid <- names(cds5p.list)[i]
    jxn <- coding_jxns[rid]
    cds.5p <- cds5p.list[[i]]
    tmp <- lapply(cds.5p, clipFivePrime, jxn)
    cds5p.truncated[[i]] <- GRangesList(tmp)
  }
  cds5p.truncated <- cds5p.truncated[ names(coding_jxns) ]
  cds5p.truncated
}

clip3prime <- function(tx, jxn){
  if(isNegStrand(tx)){
    tx2 <- use_negative_coords(tx)
    jxn["3p"] <- use_negative_coords(jxn["3p"])
  } else tx2 <- tx
  select <- start(tx) > start(jxn["3p"])
  tx[select]
}

clipThreePrime <- function(tx, jxn){
  if(isNegStrand(tx)){
    tx2 <- use_negative_coords(tx)
    jxn <- use_negative_coords(jxn)
    select <- start(tx2) > start(jxn)
    return(tx[select])
  }
  select <- start(tx) > start(jxn)
  tx[select]
}

expand.grid2 <- function(grl1, grl2){
  x <- seq_along(grl1)
  y <- seq_along(grl2)
  z <- expand.grid(x, y)
  fusion <- vector("list", nrow(z))
  nms <- rep(NA, nrow(z))
  i <- 1
  ##for(i in 1:nrow(z)){
  for(j in seq_along(grl1)){
    for(k in seq_along(grl2)){
      ##j <- z[i, 1]
      ##k <- z[i, 2]
      fivep <- grl1[[j]]
      fivep$tx_name <- rep(names(grl1)[j], length(fivep))
      threep <- grl2[[k]]
      threep$tx_name <- rep(names(grl2)[k], length(threep))
      fusion[[i]] <- c(fivep, threep)
      nms[i] <- paste(names(grl1)[j],
                      names(grl2)[k], sep="-")
      i <- i+1
    }
  }
  names(fusion) <- nms
  GRangesList(fusion)
}

arrange_clipped_tx <- function(rlist, orientation=1,
                               tx, cds){
  r <- rlist[[orientation]]
  tx.list <- rearrangementTranscripts2(linkedBins(r), tx, cds)
  if(!all(elementNROWS(tx.list)) > 0){
    return(NULL)
  }
  jxn <- seqJunction(r)
  ## all 5' transcripts
  tx5p.list <- tx.list[["5p"]]
  ## all 3' transcripts
  tx3p.list <- tx.list[["3p"]]
  tx5p.l2 <- lapply(tx5p.list, clip5prime, jxn=jxn)
  ##
  ## different transcripts can yield the exact same CDS after clipping
  ##
  tx5p.l2 <- tx5p.l2[!duplicated(tx5p.l2)]
  tx3p.l2 <- lapply(tx3p.list, clip3prime, jxn=jxn)
  tx3p.l2 <- tx3p.l2[!duplicated(tx3p.l2)]
  if(any(elementNROWS(tx5p.l2) == 0)){
    ## fiveUTRsByTranscript(txdb)
    index <- which(elementNROWS(tx5p.l2) == 0)
    for(i in index){
      nm <- names(tx5p.l2)[i]
      jxn.5p <- jxn["5p"]
      jxn.5p$cds_id <- "promoter"
      jxn.5p$cds_name <- NA
      jxn.5p$exon_rank <- seq_len(length(jxn.5p))
      ##jxn.5p$tx_name <- nm
      tmp <- GRangesList(jxn.5p)
      names(tmp) <- nm[i]
      tx5p.l2[[i]] <- tmp
    }
  }
  tx3p.l2 <- tx3p.l2[elementNROWS(tx3p.l2) > 0]
  if(length(tx3p.l2) == 0){
    browser()
  }
  fusion.grl <- expand.grid2(tx5p.l2, tx3p.l2)
  mcols(fusion.grl)$rid <- names(rlist)[orientation]
  fusion.grl
}

.fusion_GRangesList <- function(r, tx, cds){
  tx.list <- rearrangementTranscripts2(linkedBins(r), tx, cds)
  jxn <- seqJunction(r)
  ## 5' transcripts
  tx5p.list <- tx.list[["5p"]]
  ## 3' transcripts
  tx3p.list <- tx.list[["3p"]]
  tx5p.l2 <- lapply(tx5p.list, clip5prime, jxn=jxn)
  tx5p.l2 <- tx5p.l2[!duplicated(tx5p.l2)]
  tx3p.l2 <- lapply(tx3p.list, clip3prime, jxn=jxn)
  tx3p.l2 <- tx3p.l2[!duplicated(tx3p.l2)]
  fusion.grl <- expand.grid2(tx5p.l2, tx3p.l2)
  mcols(fusion.grl)$rid <- names(r)[1]
  fusion.grl
}

fusion_CDS <- function(rlist, tx, cds){
  fusion1 <- arrange_clipped_tx(rlist, orientation=1,
                                tx, cds)
  fusion2 <- arrange_clipped_tx(rlist, orientation=2,
                                tx, cds)
  if(FALSE){
    grl <- GRangesList()
    grl$orientation1 <- fusion1
    grl$orientation2 <- fusion2
  }
  ## this returned a list and not a GRangesList within this package's NAMESPACE
  fusions <- c(orientation1=GRangesList(fusion1),
               orientation2=GRangesList(fusion2))
  fusions
}

fuseCDS <- function(r, jxn, cds){
  ##tx5p.id <- gsub("promoter_", "", strsplit(jxn$tx_name, ",")[[1]])
  tx5p.id <- strsplit(jxn$tx_name, ",")[[1]]
  tx3p.id <- strsplit(jxn$"3p"$tx_name, ",")[[1]]
  cds.5p <- cds[tx5p.id]
  cds.3p <- cds[tx3p.id]
  tx5p.l2 <- lapply(cds.5p, clipFivePrime, jxn=jxn)
  tx5p <- tx5p.l2[!duplicated(tx5p.l2)]
  tx3p.l2 <- lapply(cds.3p, clipThreePrime, jxn=jxn$"3p"[1])
  tx3p <- tx3p.l2[!duplicated(tx3p.l2)]
  grl <- expand.grid2(tx5p, tx3p)
  id5p <- sapply(strsplit(names(grl), "-"), "[", 1)
  id3p <- sapply(strsplit(names(grl), "-"), "[", 2)
  tx5p <- tx5p[ id5p ]
  tx3p <- tx3p[ id3p ]
  cds.3p <- cds.3p[ id3p ]
  cds.5p <- cds.5p[ id5p ]
  names(cds.5p) <- names(cds.3p) <- names(tx5p) <- names(grl)
  list(fusion=grl,
       cds5p.tumor=GRangesList(tx5p),
       cds3p.tumor=GRangesList(tx3p),
       cds5p.ref=GRangesList(cds.5p),
       cds3p.ref=GRangesList(cds.3p))
}

#' Extract the CDS involved in each rearrangement of a RearrangementList object
#'
#' Extract for each rearrangement the full CDS of the fused sequence in the somatic genome (fusions), the partial CDS of the 5-prime (tum.5p) and 3-prime (tum.3p) transcripts that are involved in the fusion, and the full CDS of the 5-prime (ref.5p) and 3-prime transcripts in the reference genome.
#'
#' @param rlist a \code{RearrangementList}
#' @param jxns a \code{GRanges} specifying the 5-prime and 3-prime genomic regions that form a new sequence junction
#' @return a \code{List} of the CDS
#' @export
fuseCDS_Rlist <- function(rlist, jxns){
  cds <- promoterCDS(jxns, genome(jxns)[[1]])
  rlist <- rlist[jxns$rid]
  stopifnot(identical(names(rlist), names(jxns)))
  fusion.list <- vector("list", length(rlist))
  names(fusion.list) <- names(rlist)
  cds5p.list <- setNames(vector("list", length(rlist)),
                         names(fusion.list))
  tum3p.list <- setNames(vector("list", length(rlist)),
                         names(fusion.list))
  tum5p.list <- setNames(vector("list", length(rlist)),
                         names(fusion.list))
  ref5p.list <- setNames(vector("list", length(rlist)),
                         names(fusion.list))
  ref3p.list <- setNames(vector("list", length(rlist)),
                         names(fusion.list))
  for(i in seq_along(rlist)){
    tmp <- fuseCDS(rlist[[i]], jxns[i], cds)
    fusion.list[[i]] <- tmp[["fusion"]]
    tum5p.list[[i]] <- tmp[["cds5p.tumor"]]
    tum3p.list[[i]] <- tmp[["cds3p.tumor"]]
    ref5p.list[[i]] <- tmp[["cds5p.ref"]]
    ref3p.list[[i]] <- tmp[["cds3p.ref"]]
  }
  result <- list(fusions=List(fusion.list),
                 tum.5p=List(tum5p.list),
                 tum.3p=List(tum3p.list),
                 ref.5p=List(ref5p.list),
                 ref.3p=List(ref3p.list),
                 coding_junctions=jxns)
  List(result)
}

reference_protein_3prime<- function(ref_cds, genome){
  ref.seq <- DNAStringSet(lapply(getSeq(genome, ref_cds), unlist))
  translate(ref.seq, if.fuzzy.codon="solve")
}

## ENTIRE fusion protein
tumor_protein <- function(fusions, genome, all.frames=TRUE){
  tumor.seq <- extractTranscriptSeqs2(genome, fusions)
  tumor.seq.frame1 <- frame(tumor.seq, 0)
  tumor.prt.frame1 <- translate(tumor.seq.frame1, if.fuzzy.codon="solve")
  if(!all.frames) return(frame1=tumor.prt.frame1)
  tumor.seq.frame2 <- frame(tumor.seq, 1)
  tumor.seq.frame3 <- frame(tumor.seq, 2)
  tumor.prt.frame2 <- translate(tumor.seq.frame2, if.fuzzy.codon="solve")
  tumor.prt.frame3 <- translate(tumor.seq.frame3, if.fuzzy.codon="solve")
  list(frame1=tumor.prt.frame1,
       frame2=tumor.prt.frame2,
       frame3=tumor.prt.frame3)
}

#' Translate the CDS in each of the three possible reading frames
#'
#' For each fusion, we translate the CDS in each of the three possible reading frames.
#'
#' @param cds.list a \code{List} of CDS for each fusion
#' @seealso \code{\link{fuseCDS_Rlist}}
#' @export
translateCDS <- function(cds.list){
  build <- genome(cds.list[[1]][[1]])[[1]]
  bs.pkg <- paste0("BSgenome.Hsapiens.UCSC.", build)
  genome <- getBSgenome(bs.pkg)  

  fuse.list <- cds.list[["fusions"]]
  tum5p.list <- cds.list[["tum.5p"]]
  tum3p.list <- cds.list[["tum.3p"]]
  ref5p.list <- cds.list[["ref.5p"]]
  ref3p.list <- cds.list[["ref.3p"]]
  fuse.p <- lapply(fuse.list, tumor_protein, genome=genome)
  tum5p.p <- lapply(tum5p.list, tumor_protein, genome=genome,
                    all.frames=TRUE)
  tum3p.p <- lapply(tum3p.list, tumor_protein, genome=genome,
                    all.frames=TRUE)
  ref5p.p <- lapply(ref5p.list, tumor_protein, genome=genome,
                    all.frames=TRUE)
  ref3p.p <- lapply(ref3p.list, tumor_protein, genome=genome,
                    all.frames=TRUE)
  List(fusion=List(fuse.p),
       tum5p=List(tum5p.p),
       tum3p=List(tum3p.p),
       ref5p=List(ref5p.p),
       ref3p=List(ref3p.p))
}

remove_5prime_aminoAcids <- function(fusion_cds, tproteins, genome, index){
  ##tumor.cds <- fusions[[i]]
  tx.5prime <- fusion_cds$tx_name[1]
  cds.5prime <- fusion_cds[fusion_cds$tx_name==tx.5prime]
  seq.5prime <- DNAStringSet(lapply(getSeq(genome, cds.5prime), unlist))
  n.aa <- (sum(width(seq.5prime)) + 0:2)/3
  ## If the clipped 5' is tx is exactly an integer number of codons and the 3' transcript has no cds clipped, the amino acid sequence of the 3' transcript should be the same as the reference?
  select <- n.aa %% 1  == 0
  if(!any(select)) stop("did not expect to reach this point")
  drop.aa <- seq_len(n.aa [ select  ])
  tproteins[["frame1"]] <- tproteins[["frame1"]][[index]][-drop.aa]
  tproteins[["frame2"]] <- tproteins[["frame2"]][[index]][-drop.aa]
  tproteins[["frame3"]] <- tproteins[["frame3"]][[index]][-drop.aa]
  tproteins
}


is_3pTumor_sameAs_3pRef <- function(fusions, genome, cds){
  txid.3p <- unique(sapply(strsplit(names(fusions), "-"), "[", 2))
  ##ref.prt <- reference_protein_3prime(ref_cds=cds[txid.3p], genome)
  rprotein_3prime <- reference_protein_3prime(ref_cds=cds[txid.3p], genome)
  ##
  ## Assessing in-frame
  ##
  ##  1. get sequence of ENTIRE fused transcript
  ##  2. translate sequence in all 3 frames
  ##        -- remove basepairs from end of sequence so that the length of the sequence is divisible by 3 in each frame
  ##  3. Remove the amino acids in the tumor protein that correspond to amino acids from the 5' transcript.  Evaluate whether the amino acids derived from the 3' transcript of the fusion protein are a subsequence of the wild-type amino acid sequence of the reference protein
  tproteins <- tumor_protein(genome, fusions)
  L <- length(tproteins[[1]])
  ij <- expand.grid(seq_len(L), seq_along(rprotein_3prime))
  k <- 1
  in_frame <- rep(NA, nrow(ij))
  for(i in seq_len(L)){
    tumor_aa <- remove_5prime_aminoAcids(fusion_cds=fusions[[i]],
                                         tproteins=tproteins,
                                         genome=genome,
                                         index=i)
    for(j in seq_along(rprotein_3prime)){
      ref.aa <- rprotein_3prime[[j]]
      frame1.inframe <- isInFrame(query=tumor_aa[["frame1"]], ref=ref.aa)
      frame2.inframe <- isInFrame(query=tumor_aa[["frame2"]], ref=ref.aa)
      frame3.inframe <- isInFrame(query=tumor_aa[["frame3"]], ref=ref.aa)
      inframe <- c(frame1.inframe, frame2.inframe, frame3.inframe)
      in_frame[k] <- any(inframe)
      k <- k+1
    }
  }
  in_frame
}

frame <- function(dna.ss, start.pos=0){
  if(start.pos==0){
    return(divisibleBy3Seq(dna.ss))
  }
  dna.list <- lapply(dna.ss, function(x, i){
    exclude <- seq_len(i)
    x[-exclude]
  }, i = start.pos)
  divisibleBy3Seq(DNAStringSet(dna.list))
}


candidateFusions <- function(rlist, cds){
  lb <- linkedBins(rlist)
  region1 <- lb
  region2 <- linkedTo(lb)
  cds.region1 <- overlapsAny(region1, cds, ignore.strand=TRUE, maxgap=5000)
  ## the 3' sequence junction must be in the transcript
  cds.region2 <- overlapsAny(region2, cds, ignore.strand=TRUE, maxgap=0)
  cds.region1 & cds.region2
}




list_tx_by_query <- function(query, subject, ignore.strand=TRUE, promoter=TRUE, ...){
  hits <- findOverlaps(query, subject, ignore.strand=TRUE, ...)
  if(length(hits) == 0) return(NULL)
  if(promoter){
    tx.names <- paste0("promoter_", subject$tx_name[subjectHits(hits)])
  } else {
    tx.names <- subject$tx_name[subjectHits(hits)]
  }
  tx.list <- split(tx.names, queryHits(hits))
  tx.list <- lapply(tx.list, paste, collapse=",")
  tx <- unlist(tx.list)
  tx
}

annotate5prime <- function(jxn.5p, txdb){
  tx <- transcripts(txdb)
  tx <- tx[grep("NM_", tx$tx_name)]
  prm <- promoters(txdb, upstream=5000)
  prm <- prm[grep("NM_", prm$tx_name)]
  prm <- list_tx_by_query(query=jxn.5p, subject=prm, promoter=TRUE)
  prm.index <- as.integer(names(prm))
  jxn.5p$tx_name <- ""
  jxn.5p$tx_name[prm.index] <- prm
  ##
  ##
  promoter.index <- which(jxn.5p$tx_name != "")
  tx.nms <- list_tx_by_query(query=jxn.5p, subject=tx, promoter=FALSE)
  tx.nms <- tx.nms[ ! names(tx.nms) %in% as.character(promoter.index) ]
  nonpromoter.index <- as.integer(names(tx.nms))
  jxn.5p$tx_name[nonpromoter.index] <- tx.nms
  jxn.5p
}

annotate3prime <- function(jxn.3p, txdb){
  tx <- transcripts(txdb)
  tx <- tx[grep("NM_", tx$tx_name)]
  tx.nms <- list_tx_by_query(query=jxn.3p, subject=tx, promoter=FALSE)
  nonpromoter.index <- as.integer(names(tx.nms))
  jxn.3p$tx_name <- ""
  jxn.3p$tx_name[nonpromoter.index] <- tx.nms
  jxn.3p
}

#' Select sequence junctions in which both the five-prime and three-prime genomic regions occur within a transcript
#'
#' For the 3-prime genomic region, the sequence junction must occur within a transcript (usually intronic).  For the 5-prime sequence junction, the genomic region can either be within 5-kb upstream of the transcript (defined as the promoter) or within the transcript.
#'
#' @param jxns a \code{GRanges} representation of the sequence junctions
#' @param build a character-string
#' @export
#' @examples
#' # TODO
codingJunctions <- function(jxns, build){
  txdb <- loadTxdbTranscripts(build, seqlevels(jxns))[["txdb"]]
  fivep <- annotate5prime(jxns, txdb)
  threep <- annotate3prime(jxns$"3p", txdb)
  fivep$"3p" <- threep
  names(fivep) <- fivep$rid
  select <- fivep$tx_name != "" & threep$tx_name != ""
  jxns <- fivep[ select ]
  genome(jxns) <- build
  jxns
}

referenceDna3p <- function(jxn, cds){
  txid.3p <- strsplit(jxn$"3p"$tx_name, ",")[[1]]
  ref.tx <- cds[txid.3p]
  ref.seq <- DNAStringSet(lapply(getSeq(genome, ref.tx), unlist))
  ref.seq
}

referenceDna5p <- function(jxn, cds){
  txid.5p <- strsplit(jxn$tx_name, ",")[[1]]
  ref.tx <- cds[txid.5p]
  ref.seq <- DNAStringSet(lapply(getSeq(genome, ref.tx), unlist))
  ref.seq
}

## number of amino acids from 5-prime transcript
number5pAminoAcids <- function(cds.5prime, genome){
  n.codons <- rep(NA, length(cds.5prime))
  for(i in seq_along(cds.5prime)){
    seq.5prime <- DNAStringSet(unlist(DNAStringSet(lapply(getSeq(genome, cds.5prime[[i]]), unlist))))
    n.tmp <- (width(seq.5prime) + 0:2)
    n.seq <- n.tmp[n.tmp %% 3 == 0]
    n.codons[i] <- n.seq/3
  }
  n.codons
}

number5pAminoAcids_list <- function(cds.list, genome){
  number_aa <- vector("list", length(cds.list))
  names(number_aa) <- names(cds.list)
  for(i in seq_along(cds.list)){
    cds.grl <- cds.list[[i]]
    nr <- elementNROWS(cds.grl)
    tmp <- rep(NA, length(cds.grl))
    tmp[ nr == 0 ] <- 0
    if(any(nr > 0)){
      tmp [ nr > 0 ] <- number5pAminoAcids(cds.grl[ nr > 0 ], genome)
    }
    number_aa[[i]] <- tmp
  }
  number_aa
}

aaThreePrime <- function(tprotein, fiveprimeCount){
  if(all(fiveprimeCount == 0)) return(tprotein)
  for(i in 1:3){
    prt.list <- tprotein[[i]]
    for(j in seq_along(prt.list)){
      n.AA <- fiveprimeCount[j]
      if(n.AA == 0) next()
      drop.aa <- seq_len(n.AA)
      prt.list[[j]] <- prt.list[[j]][-drop.aa]
    }
    tprotein[[i]] <- prt.list
  }
  tprotein
}

aaThreePrime_list <- function(tprotein.list, number.aa.list){
  stopifnot(identical(names(tprotein.list), names(number.aa.list)))
  for(i in seq_along(tprotein.list)){
    n.aa <- number.aa.list[[i]]
    tprotein <- tprotein.list[[i]]
    tprotein.list[[i]] <- aaThreePrime(tprotein, n.aa)
  }
  tprotein.list
}

aaFivePrime <- function(tprotein, fiveprimeCount){
  select.aa <- seq_len(fiveprimeCount-1)
  for(i in 1:3){
    prt.list <- tprotein[[i]]
    for(j in seq_along(prt.list)){
      prt.list[[j]] <- prt.list[[j]][select.aa]
    }
    tprotein[[i]] <- prt.list
  }
  tprotein
}

#' Create the full CDS of transcripts for those in which only the promoter of the transcript is involved in a fusion
#'
#' Rearrangements in which the 5-prime genomic region involves only the promoter are named promoter_TX_NAME, where TX_NAME is the name of a transcript. This function returns a GRangesList of CDS, including CDS with name promoter_TX_NAME.
#'
#' @param coding_jxns a \code{GRanges} representation of the sequence junctions that occur in transcripts, possibly intronic.
#' @param cds a \code{GRangesList} of all CDS created by \code{cdsBy}
#' @seealso \code{\link{codingJunctions}} \code{\link[GenomicFeatures]{cdsBy}}
#' @return \code{GRangesList}
promoterCDS <- function(coding_jxns, build){
  tx.cds <- loadTxdbTranscripts(build, seqlevels(coding_jxns))
  cds <- tx.cds[["cds"]]
  tx.names <- coding_jxns$tx_name
  index <- grep("promoter", tx.names)
  if(length(index) < 1) return(cds)
  prm.names <- tx.names[index]
  prm.names <- unlist(strsplit(prm.names, ","))
  promoter.tx <- gsub("promoter_", "", prm.names)
  cds.promoter <- cds[promoter.tx]
  names(cds.promoter) <- prm.names
  c(cds, cds.promoter)
}

fivePrimeTranscriptNames <- function(jxns){
  tx.list <- strsplit(jxns$tx_name, ",")
  tx.nms <- unlist(tx.list)
  rid <- rep(names(jxns), elementNROWS(tx.list))
  names(tx.nms) <- rid
  tx.nms
}

## Assumes fusion and tumor.5p are in the same frame
partitionAASequence_oneframe <- function(fusion, tumor.5p){
  ##
  ## number of amino acids derived from 5' sequence
  ##
  nAA <- lapply(tumor.5p, width)
  fusion.5p <- fusion.3p <- fusion
  for(j in seq_along(fusion)){
    Nvec <- nAA[[j]]
    Pvec <- Pvec.5p <- Pvec.3p <- fusion[[j]]
    for(k in seq_along(Pvec)){
      N <- Nvec[k]
      P <- Pvec[[k]]
      L <- length(P)
      Pvec.3p[[k]] <- subseq(P, start=N+1, end=L)
      ## when N=0, this will be a 0-letter AAString
      Pvec.5p[[k]] <- subseq(P, start=1, end=N)
    }
    fusion.5p[[j]] <- Pvec.5p
    fusion.3p[[j]] <- Pvec.3p
  }
  list("3p"=List(fusion.3p),
       "5p"=List(fusion.5p),
       fusion=List(fusion))
}

#' Partition the amino acid sequence of the fusion protein into the sequences derived from the 5-prime and 3-prime transcripts
#'
#' By partioning the amino acid sequence of the fusion into its 5-prime and 3-prime parts, we can compare the 5-prime partition to the amino acid sequence of the 5-prime reference protein and the 3-prime partition to the amino acid sequence of the 3-prime reference protein.
#'
#' @seealso \code{\link{translateCDS}}
#' @param proteins a \code{List} of amino acid sequences
#' @export
partitionAASequence <- function(proteins){
  fusion.frame1 <- lapply(proteins[["fusion"]], "[[", 1)
  tumor.5p.frame1 <- lapply(proteins[["tum5p"]], "[[", 1)
  frame1.part <- partitionAASequence_oneframe(fusion=fusion.frame1,
                                     tumor.5p=tumor.5p.frame1)
  fusion.frame2 <- lapply(proteins[["fusion"]], "[[", 2)
  tumor.5p.frame2 <- lapply(proteins[["tum5p"]], "[[", 2)
  frame2.part <- partitionAASequence_oneframe(fusion=fusion.frame2,
                                     tumor.5p=tumor.5p.frame2)
  fusion.frame3 <- lapply(proteins[["fusion"]], "[[", 3)
  tumor.5p.frame3 <- lapply(proteins[["tum5p"]], "[[", 3)
  frame3.part <- partitionAASequence_oneframe(fusion=fusion.frame3,
                                              tumor.5p=tumor.5p.frame3)
  List(frame1=List(frame1.part),
       frame2=List(frame2.part),
       frame3=List(frame3.part))
}

positionStopCodons <- function(aa){
  ## stop codon at the end is expected, so remove the last AA 
  aa <- subseq(aa, 1, length(aa)-1)
  pos <- gregexpr("\\*", as.character(aa))[[1]]
  as.numeric(pos)
}

stopPositions_oneframe<- function(fusion){
  stop.positions <- setNames(vector("list", length(fusion$"5p")),
                             names(fusion$"5p"))
  tumor.protein <- fusion[["fusion"]]
  for(i in seq_along(tumor.protein)){
    Tvec <- tumor.protein[[i]]
    ## loop over proteins from different transcripts for a particular fusion
    stops <- rep(NA, length(Tvec))
    for(j in seq_along(Tvec)){
      tmp <- positionStopCodons(Tvec[[j]])
      stops[j] <- paste(tmp, collapse=",")
    }
    stop.positions[[i]] <- stops
  }
  List(stop.positions)
}

#' Assess whether there are any premature stop codons 
#'
#' For a fusion to be in-frame, we require that the only stop codon is at the 3-prime end of the fusion amino acid sequence.
#'
#' @param fusion  a \code{List} created by \code{partitionAASequence}
#' @seealso \code{\link{partitionAASequence}}
#' @return a \code{LogicalList}
#' @export
noPrematureStop <- function(fusion){
  nostop.list <- setNames(vector("list", 3),
                          paste0("frame", 1:3))
  for(i in 1:3){
    tmp <- stopPositions_oneframe(fusion[[i]])
    nostop.list[[i]] <- tmp == "-1"
  }
  nostop.list
}

#' Combine the in-frame and no-premature-stop LogicalLists
#'
#' This function combines two \code{LogicalList}s for each reading frame, creating a new \code{LogicalList} list that is \code{TRUE} if the rearrangement is in-frame and has no premature stop codons.
#'
#' @param nostop.list a length=3 \code{List} of \code{LogicalList} elements
#' @param inframe.list a length=3 \code{List} of \code{LogicalList} elements
#' @seealso \code{\link{inFrameList}} \code{\link{noPrematureStop}}
#' @export
inFrameNoStop <- function(nostop.list, inframe.list){
  inframe_nostop <- setNames(vector("list", 3), names(nostop.list))
  frame.names <- paste0("frame", 1:3)
  nostop.list <- nostop.list[frame.names]
  inframe.list <- inframe.list[frame.names]
  for(i in 1:3){
    inframe_nostop[[i]] <- nostop.list[[i]] & inframe.list[[i]]
  }
  inframe_nostop
}

isInFrame_oneframe <- function(fusion, reference){
  inframe <- setNames(vector("list", length(fusion$"5p")),
                      names(fusion$"5p"))
  fusion.5p <- fusion$"5p"
  fusion.3p <- fusion$"3p"
  ref.5p <- reference$"5p"
  ref.3p <- reference$"3p"
  for(i in seq_along(fusion.5p)){
    T5p <- fusion.5p[[i]]
    T3p <- fusion.3p[[i]]
    R5p <- ref.5p[[i]]
    R3p <- ref.3p[[i]]
    inFr <- rep(NA, length(T5p))
    for(j in seq_along(T5p)){
      ##stops.5p <- positionStopCodons(T5p[[j]])
      ##stops.3p <- positionStopCodons(T3p[[j]], fiveprime=FALSE)
      inFr[j] <- isInFrame(T5p[[j]], R5p[[j]]) &&
        isInFrame(T3p[[j]], R3p[[j]])
    }
    inframe[[i]] <- inFr
  }
  LogicalList(inframe)
}


#' Evaluate whether fusion amino acid sequence is a subsequence of the reference
#' 
#' For each fusion and for each possible read frame, we check whether the amino acid sequence derived from the 5- and 3-prime transcripts of the fusion are subsequences of the full amino acid sequence of the reference genome translated in the same frame
#'
#' @param fusion.frames a \code{List} created by \code{partitionAASequence}
#' @param ref.frames a \code{List} created by \code{organizeReferenceByFrame}
#' @seealso \code{\link{partitionAASequence}} \code{\link{organizeReferenceByFrame}}
#' @export
inFrameList <- function(fusion.frames, ref.frames){
  inframe <- vector("list", 3)
  for(i in 1:3){
    inframe[[i]] <- isInFrame_oneframe(fusion.frames[[i]], ref.frames[[i]])
  }
  ##any_inframe <- inframe[[1]] | inframe[[2]] | inframe[[3]]
  ##any_inframe
  names(inframe) <- paste0("frame", 1:3)
  inframe
}

#' Organize reference amino acid sequences by frame
#'
#' @param proteins an object created by \code{translateCDS}
#' @return a \code{List} of the reference amino acid sequences for each of the three possible reading frames
#' @export
#' @seealso \code{\link{translateCDS}}
organizeReferenceByFrame <- function(proteins){
  ref5p.frame1 <- lapply(proteins[["ref5p"]], "[[", 1)
  ref5p.frame2 <- lapply(proteins[["ref5p"]], "[[", 2)
  ref5p.frame3 <- lapply(proteins[["ref5p"]], "[[", 3)

  ref3p.frame1 <- lapply(proteins[["ref3p"]], "[[", 1)
  ref3p.frame2 <- lapply(proteins[["ref3p"]], "[[", 2)
  ref3p.frame3 <- lapply(proteins[["ref3p"]], "[[", 3)

  ref.frame1 <- list("5p"=ref5p.frame1, "3p"=ref3p.frame1)
  ref.frame2 <- list("5p"=ref5p.frame2, "3p"=ref3p.frame2)
  ref.frame3 <- list("5p"=ref5p.frame3, "3p"=ref3p.frame3)
  List(frame1=List(ref.frame1),
       frame2=List(ref.frame2),
       frame3=List(ref.frame3))
}

#' Create a table of in-frame fusions
#'
#' @param fusions a \code{List} created by \code{validFusions}
#' @param jxns a \code{GRanges} object of the sequence junctions indicating the 5-prime and 3-prime genomic regions
#' @seealso \code{\link{validFusions}}
#' @export
fusionTable2 <- function(fusions){
  jxns <- fusions[["coding_junctions"]]
  build <- genome(jxns)[[1]]
  tx <- loadTx(build)
  genes <- setNames(tx$gene_name, tx$tx_name)

  cds.5p <- fusions[["cds.5p"]]
  cds.3p <- fusions[["cds.3p"]]
  ref.5p <- fusions[["ref.5p"]]
  ref.3p <- fusions[["ref.3p"]]

  tx.5prime <- unlist(lapply(cds.5p, names))
  tx.5prime <- sapply(strsplit(tx.5prime, "-"), "[", 1)
  tx.3prime <- unlist(lapply(cds.3p, names))

  nexons <- List(ref.3p=unlist(lapply(ref.3p, elementNROWS)),
                 tum.3p=unlist(lapply(cds.3p, elementNROWS)),
                 tum.5p=unlist(lapply(cds.5p, elementNROWS)))
  first.5p <- pmin(1, nexons$tum.5p)
  last.5p <- nexons$tum.5p
  exons.5p <- paste0("exons: ", first.5p, "-", last.5p)

  last.3p <- nexons[["ref.3p"]]
  first.3p <- last.3p - nexons[["tum.3p"]] + 1
  exons.3p <- paste0("exons: ", first.3p, "-", last.3p)
  exons.3p[first.3p > last.3p] <- "none"

  nms <- rep(names(cds.5p), elementNROWS(cds.5p))
  jxns <- jxns[ nms ]
  fivep.jxn <- paste0(chromosome(jxns), ":", start(jxns))
  threep.jxn <- paste0(chromosome(jxns$"3p"), ":", start(jxns$"3p"))
  fivep.strand <- as.character(strand(jxns))
  threep.strand <- as.character(strand(jxns$"3p"))
  tx.nms.5p <- gsub("promoter_", "", tx.5prime)
  genes.5p <- genes[tx.nms.5p]
  genes.3p <- genes[tx.3prime]
  aa.5p <- fusions$tum_aa.5p
  aa.5p.end <- unlist(lapply(aa.5p, lengths))
  aa.5p.start <- pmin(aa.5p.end, 1)
  ##
  ## To determine the starting index of the 3prime gene, we need to know the length of the reference amino acid
  ##
  ref.3p <- fusions$ref_aa.3p
  tum.3p.end <- unlist(lapply(ref.3p, lengths))
  tum.3p <- fusions$tum_aa.3p
  tum.3p.length <- unlist(lapply(tum.3p, lengths))
  tum.3p.start <- tum.3p.end - tum.3p.length + 1
  df <- DataFrame(gene.5prime=genes.5p,
                  gene.3prime=genes.3p,
                  tx.5prime=tx.5prime,
                  tx.3prime=tx.3prime,
                  junction.5prime=fivep.jxn,
                  junction.3prime=threep.jxn,
                  exons.5prime=exons.5p,
                  exons.3prime=exons.3p,
                  strand.5prime=fivep.strand,
                  strand.3prime=threep.strand,
                  rearrangement.id=names(jxns),
                  aa.5prime.start=aa.5p.start,
                  aa.5prime.end=aa.5p.end,
                  aa.3prime.start=tum.3p.start,
                  aa.3prime.end=tum.3p.end)
  df
}

#' Create a List of in-frame fusions with no premature stops
#'
#' @param fusion.list a \code{List} created by \code{partitionAASequence}
#' @param cds.list a \code{List} created by \code{fuseCDS_Rlist}
#' @param inframe.nostop a \code{List} of \code{LogicalList}s.  True values indicate fusions that are both in-frame and have no premature stop codons.
#' @param ref.frames list of amino acid squences for reference proteins, organized by frame
#' @seealso \code{\link{partitionAASequence}} \code{\link{fuseCDS_Rlist}}
#' @export
validFusions <- function(fusion.list, cds.list, inframe.nostop, ref.frames){
  nms <- paste0("frame", 1:3)
  fusion.list2 <- setNames(vector("list", 3), nms)
  ref.3p <- ref.5p <- aa.3p <- aa.5p <- setNames(vector("list", 3), nms)
  fusion.list <- fusion.list[nms]
  inframe.list <- inframe.nostop[nms]
  for(i in 1:3){
    fusions <- fusion.list[[i]][["fusion"]]
    is.valid <- inframe.nostop[[i]]
    fusions <- fusions[ is.valid ]
    fusions <- fusions[ elementNROWS(fusions) > 0 ]
    fusion.list2[[i]] <- fusions

    temp.5p <- as(ref.frames[[i]][["5p"]], "List")
    temp.5p <- temp.5p[names(is.valid)]
    temp.3p <- as(ref.frames[[i]][["3p"]], "List")
    temp.3p <- temp.3p[names(is.valid)]
    temp.5p <- temp.5p[ is.valid ]
    temp.3p <- temp.3p[ is.valid ]
    ref.5p[[i]] <- temp.5p[ elementNROWS(temp.5p) > 0 ]
    ref.3p[[i]] <- temp.3p[ elementNROWS(temp.3p) > 0 ]

    tmp.5p <- fusion.list[[i]][["5p"]]
    tmp.3p <- fusion.list[[i]][["3p"]]
    tmp.5p <- tmp.5p [ is.valid ]
    tmp.3p <- tmp.3p [ is.valid ]
    aa.5p[[i]] <- tmp.5p[ elementNROWS(tmp.5p) > 0 ]
    aa.3p[[i]] <- tmp.3p[ elementNROWS(tmp.3p) > 0 ]
  }
  nms.frame1 <- names(fusion.list2[["frame1"]])
  nms.frame2 <- names(fusion.list2[["frame2"]])
  nms.frame3 <- names(fusion.list2[["frame3"]])
  nms.frame2 <- nms.frame2[!nms.frame2 %in% nms.frame1]
  nms.frame3 <- nms.frame3[!nms.frame3 %in% c(nms.frame1, nms.frame2)]
  fusion.aa <- c(fusion.list2[["frame1"]],
                 fusion.list2[["frame2"]][nms.frame2],
                 fusion.list2[["frame3"]][nms.frame3])
  aa.5p <- c(aa.5p[["frame1"]],
             aa.5p[["frame2"]][nms.frame2],
             aa.5p[["frame3"]][nms.frame3])
  aa.3p <- c(aa.3p[["frame1"]],
             aa.3p[["frame2"]][nms.frame2],
             aa.3p[["frame3"]][nms.frame3])
  ref_aa.5p <- c(ref.5p[["frame1"]],
                 ref.5p[["frame2"]][nms.frame2],
                 ref.5p[["frame3"]][nms.frame3])
  ref_aa.3p <- c(ref.3p[["frame1"]],
                 ref.3p[["frame2"]][nms.frame2],
                 ref.3p[["frame3"]][nms.frame3])
  fusion.cds <- cds.list[["fusions"]]
  fusion.cds <- fusion.cds[ names(fusion.aa) ]
  is.valid <- inframe.nostop[[1]]
  cds.5p <- cds.list[["tum.5p"]]
  cds.5p <- cds.5p[ names(is.valid) ]
  cds.5p <- cds.5p[ is.valid ]
  cds.5p <- cds.5p[ elementNROWS(cds.5p) > 0 ]

  cds.3p <- cds.list[["tum.3p"]]
  cds.3p <- cds.3p[ names(is.valid) ]
  cds.3p <- cds.3p[ is.valid ]
  cds.3p <- cds.3p[ elementNROWS(cds.3p) > 0 ]

  ref.5p <- cds.list[["ref.5p"]][ is.valid ]
  ref.5p <- ref.5p[ elementNROWS(ref.5p) > 0 ]
  ref.3p <- cds.list[["ref.3p"]][ is.valid ]
  ref.3p <- ref.3p[ elementNROWS(ref.3p) > 0 ]
  if(!all(names(aa.5p) %in% names(cds.5p))){
    nm <- names(aa.5p)[! names(aa.5p) %in% names(cds.5p)]
    tmp <- cds.list[["tum.5p"]][nm]
    is.valid <- inframe.nostop[[1]]
    is.valid <- is.valid[[nm]]
    if(!any(is.valid)){
      is.valid <- inframe.nostop[[2]]
      is.valid <- is.valid[[nm]]
    }
    if(!any(is.valid)){
      is.valid <- inframe.nostop[[3]]
      is.valid <- is.valid[[nm]]
    }
    tmp[[1]] <- tmp[[1]][ is.valid ]
    cds.5p <- c(cds.5p, tmp)
    cds.5p <- cds.5p[ names(aa.5p) ]

    tmp <- cds.list[["tum.3p"]][nm]
    tmp[[1]] <- tmp[[1]][ is.valid ]
    cds.3p <- c(cds.3p, tmp)
    cds.3p <- cds.3p[ names(aa.3p) ]

    tmp <- cds.list[["ref.5p"]][ nm ]
    tmp[[1]] <- tmp[[1]][ is.valid ]
    ref.5p <- c(ref.5p, tmp)
    ref.5p <- ref.5p[ names(aa.5p) ]

    tmp <- cds.list[["ref.3p"]][ nm ]
    tmp[[1]] <- tmp[[1]][ is.valid ]
    ref.3p <- c(ref.3p, tmp)
    ref.3p <- ref.3p[ names(aa.3p) ]
  }
  stopifnot(identical(lengths(cds.5p), lengths(aa.5p)))
  stopifnot(identical(lengths(cds.3p), lengths(aa.3p)))
  stopifnot(identical(lengths(ref.5p), lengths(cds.5p)))
  stopifnot(identical(lengths(ref.3p), lengths(cds.3p)))
  List(amino.acids=fusion.aa,
       tum_aa.5p=aa.5p,
       tum_aa.3p=aa.3p,
       ref_aa.5p=ref_aa.5p,
       ref_aa.3p=ref_aa.3p,
       cds=fusion.cds,
       cds.5p=cds.5p,
       cds.3p=cds.3p,
       ref.5p=ref.5p,
       ref.3p=ref.3p,
       coding_junctions=cds.list[["coding_junctions"]])
}

## #' Puts txdb, BSgenome, transcripts, and CDS in global environment
## #'
## #' @param build genome build -- only 'hg19' and 'hg18' supported
## #' @export
## #' @examples
## #' rm(list=ls())
## #' loadGenomeData()
## #' exists(c("cds", "txdb", "tx", "genome"))
## loadGenomeData <- function(build="hg19"){
##   tx.cds <- loadTxdbTranscripts(build)
##   tx <- tx.cds[["transcripts"]]
##   cds <- tx.cds[["cds"]]
##   txdb <- tx.cds[["txdb"]]
## 
##   orgdb <- loadOrgDb()
##   bs.pkg <- paste0("BSgenome.Hsapiens.UCSC.", build)
##   genome <- getBSgenome(bs.pkg)
## 
##   assign("cds", cds, envir=globalenv())
##   assign("genome", genome, envir=globalenv())
##   assign("txdb", txdb, envir=globalenv())
##   assign("tx", tx, envir=globalenv())
## }


#' Convert tumor amino acid ranges to a GRanges object
#'
#' Convert tabled amino acid ranges to a GRanges object to facilitate the overlap with protein domains from Uniprot database.
#' @param tab a \code{DataFrame} of fusions as obtained from \code{fusionTable2}
#' @examples
#'  extdata <- system.file("extdata", package="trellis")
#'  fusions <- readRDS(file.path(extdata, "valid_fusions.rds"))
#'  tab <- fusionTable2(fusions)
#'  up <- readRDS(file.path(extdata, "uniprot.rds"))
#'  up2 <- uniprotFeatures(up, tab, strwrap.width=20)
#'  ## check that none of the short descriptions are too long
#'  nchar.desc <- nchar(up2$short.desc)
#'  gene1.in.up <- tab$gene.5prime %in% up2$hugo
#'  gene2.in.up <- tab$gene.3prime %in% up2$hugo
#'  bothin.up <- gene1.in.up & gene2.in.up
#'  tab <- tab[bothin.up, ]
#'  ## the aa_len field of the amino acid in uniprot might have a zero-based index
#'  ## -- it is one smaller than the length we've recorded
#'  ## coerce to GRanges
#'  tumor_aa_ranges <- aa_granges(tab)
#'  identical(as.character(seqnames(tumor_aa_ranges)),
#'            rep(c("ERBB4", "IKZF2"), each=2))
#'  domain_aa_ranges <- GRanges(up2$hugo, IRanges(up2$start, up2$end),
#'                              chromosome=up2$seqnames,
#'                              description=up2$description,
#'                              short.desc=up2$short.desc,
#'                              aa_len=up2$aa_len)
#'  domains <- subsetByOverlaps(domain_aa_ranges, tumor_aa_ranges, type="within")
#'  domains
#' @export
aa_granges <- function(tab){
   jxn.5prime <- strsplit(tab$junction.5prime, ":")
   jxn.3prime <- strsplit(tab$junction.3prime, ":")
   chr.5prime <- sapply(jxn.5prime, "[", 1)
   chr.3prime <- sapply(jxn.3prime, "[", 1)

   starts.5prime <- tab$aa.5prime.start
   ends.5prime <- tab$aa.5prime.end
   starts.3prime <- tab$aa.3prime.start
   ends.3prime <- tab$aa.3prime.end
   GRanges(c(tab$gene.5prime, tab$gene.3prime), IRanges(c(starts.5prime, starts.3prime),
                                                        c(ends.5prime, ends.3prime)),
           transcript=c(tab$tx.5prime, tab$tx.3prime),
           chromosome=c(chr.5prime, chr.3prime))
}
