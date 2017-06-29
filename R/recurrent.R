#' Annotate table of amplified genes with linked drivers and the number of
#' samples for which a driver is linked
#'
#' @param rcnvs a table of recurrent CNVs
#' @param grl a \code{GRangesList} of the segmented copy number (each element is a sample)
#'
#' @seealso See \code{\link{recurrentDrivers}} for making a table of recurrent amplicons
#' @export
annotateRecurrent <- function(rcnvs, grl){
  df_ext <- data.frame(gene=rep(rcnvs$gene, rcnvs$freq),
                       chr=rep(rcnvs$chr, rcnvs$freq),
                       start=rep(rcnvs$start, rcnvs$freq),
                       end=rep(rcnvs$end, rcnvs$freq),
                       id=unlist(strsplit(as.character(rcnvs$id), ",")))
  df_ext$driver <- NA
  names(grl) <- gsub(".bam", "", names(grl))
  dflist <- split(df_ext, df_ext$id)
  grl <- grl[names(dflist)]

  for(k in seq_along(dflist)){
    amps <- grl[[k]]
    dat <- dflist[[k]]
    dat_gr <- GRanges(dat$chr, IRanges(dat$start, dat$end))
    hits <- findOverlaps(dat_gr, amps)
    i <- queryHits(hits)
    j <- subjectHits(hits)
    j <- j[!duplicated(i)]
    i <- i[!duplicated(i)]
    dat$log2_ratio <- NA
    dat$log2_ratio[i] <- amps$seg.mean[j]
    dat$linked_driver <- NA
    dat$linked_driver[i] <- amps$driver[j]
    dflist[[k]] <- dat
  }
  df_ext2 <- do.call(rbind, dflist)
  ld <- df_ext2$linked_driver
  ld[ ld == "" ] <- NA
  df_ext2$linked_driver <- ld
  df_ext2 <- df_ext2[order(df_ext2$gene), ]
  rownames(df_ext2) <- NULL

  dflist2 <- split(df_ext2, df_ext2$gene)
  df.gene <- vector("list", length(df_ext2))
  for(i in seq_along(dflist2)){
    dat <- dflist2[[i]]
    x <- as.character(dat$linked_driver)
    x <- unique(unlist(strsplit(x, ", ")))
    if(!all(is.na(x))){
      x <- x[ !is.na(x) ]
      xx <- paste(unique(unlist(x)), collapse=",")
    } else xx <- NA
    nlinked <- sum(!is.na(dat$linked_driver))
    dat2 <- data.frame(gene=names(dflist2)[i],
                       log2_ratio=mean(dat$log2_ratio, na.rm=TRUE),
                       id=paste(dat$id, collapse=","),
                       linked_drivers=xx,
                       nlinked=nlinked,
                       freq=nrow(dat))
    df.gene[[i]] <- dat2
  }
  df.gene <- do.call(rbind, df.gene)
  df.gene <- df.gene[order(df.gene$freq, decreasing=TRUE), ]
  df.gene
}

labelGeneName <- function(tx, gr, maxgap){
  hits <- findOverlaps(tx, gr, maxgap=maxgap)
  gene_name <- rep(NA, length(gr))
  gene_name[subjectHits(hits)] <- tx$gene_name[queryHits(hits)]
  gene_name
}

## #' Summarize frequency of a sample in a GRangesList
## #'
## #' Makes a gene-level table with frequency that gene appears and a vector of
## #' comma-delimited sample identifiers with the gene. Useful for summarizing
## #' recurrence of deletions or amplicons by gene.
## #'
## #' @param grl a \code{GRangesList}.  Assumes that each element of the \code{GRangesList} has a \code{id} and a \code{gene_name} field.
## #' @seealso \code{\link{recurrentDeletions}}
## #' @examples
## #'   genes <- GRanges("chr1", IRanges(5, 10))
## #'   genes$gene_name <- "a"
## #'   gr1 <- GRanges(rep("chr1", 2), IRanges(c(4, 8), c(6, 10)), id=rep("id1", 2))
## #'   gr2 <- GRanges("chr1", IRanges(3, 9), id="id2")
## #'   grl <- GRangesList(id1=gr1, id2=gr2)
## #'   grl <- recurrentDeletions(genes, grl, maxgap=5)
## #'   tab <- summarizeGeneFreq(grl)
summarizeGeneFreq <- function(grl){
  g <- unlist(grl)
  genes <- unique(g$gene_name)
  id.list <- lapply(split(g$id, g$gene_name), unique)
  ids <- sapply(id.list, paste, collapse=",")
  freq <- elementNROWS(id.list)
  tab <- data.frame(gene=names(freq),
                    frequency=freq,
                    ids=ids)
}

#' Aggregate deletions to the gene-level when evaluating recurrence
#'
#'
#' @param tx a \code{GRanges} object of transcripts
#' @param grl a \code{GRangesList} of deletions -- each element is a sample
#' @param maxgap length-one integer vector passed to \code{findOverlaps}
#' @export
#' @examples
#'   genes <- GRanges("chr1", IRanges(5, 10))
#'   genes$gene_name <- "a"
#'   gr1 <- GRanges(rep("chr1", 2), IRanges(c(4, 8), c(6, 10)), id=rep("id1", 2))
#'   gr2 <- GRanges("chr1", IRanges(3, 9), id="id2")
#'   grl <- GRangesList(id1=gr1, id2=gr2)
#'   recurrentDeletions(genes, grl, maxgap=5)
recurrentDeletions <- function(tx, grl, maxgap=5e3){
  gr <- unlist(grl)
  gr$gene_name <- labelGeneName(tx, gr, maxgap)
  gene.freq <- table(gr$gene_name)
  recurrent.genes <- names(gene.freq)[gene.freq >= 2]
  gr <- gr[gr$gene_name %in% recurrent.genes]
  grl <- split(gr, gr$id)
  grl
}
