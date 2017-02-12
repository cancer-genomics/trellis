#' Annotate table of amplified genes with linked drivers and the number of
#' samples for which a driver is linked
#'
#' @param ramps a table of recurrent amplicons
#' @param amp_grl a \code{GRangesList} of amplicons (each element is a sample)
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

#' Aggregate deletions to the gene-level when evaluating recurrence
#'
#' @param tx a \code{GRanges} object of transcripts
#' @param grl a \code{GRangesList} of deletions -- each element is a sample
#' @param maxgap length-one integer vector passed to \code{findOverlaps}
#' @export
recurrentDeletions <- function(tx, grl, maxgap=5e3){
  browser()
  gr <- unlist(grl)
  is_overlap <- is_overlap[ cnts > 1, ] 
  cnts <- cnts[ cnts > 1 ]
  ids <- apply(is_overlap, 1, function(is_amplicon, id) {
    paste(id[is_amplicon], collapse=",")
  }, id=gsub(".bam", "", colnames(is_overlap)))
  result <- data.frame(gene = tx$gene_name, freq=as.integer(cnts), id=ids)
  ##
  ## Gene coordinates
  ##
  tx2 <- tx[tx$gene_name %in% result$gene]
  tx2.list <- GRangesList(sapply(split(tx2, tx2$gene_name), reduce))
  tx2 <- unlist(tx2.list)
  tx2 <- tx2[result$gene]
  stopifnot(identical(names(tx2), as.character(result$gene)))
  result$chr <- chromosome(tx2)
  result$start <- start(tx2)
  result$end <- end(tx2)
  result <- result[, c("gene", "chr", "start", "end", "freq", "id")]
  rownames(result) <- NULL
  result[order(result$freq, decreasing=TRUE), ]
}
