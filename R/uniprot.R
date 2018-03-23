reduceProtein <- function(gup, hugo, strwrap.width){
  gup <- gup[gup$hugo==hugo]
  short.desc <- as.character(gup$shrt_desc)
  gup$short.desc <- short.desc
  red <- reduce(gup, with.revmap=TRUE)
  revmap <- red$revmap
  gup2 <- relist(gup[unlist(revmap)], revmap)
  ## if short description is null, remove
  gup3 <- endoapply(gup2, function(g){
    g2 <- g[g$short.desc != "NULL"]
    if(length(g2) == 0){
      g2 <- g
    }
    g2
  })
  gup3 <- unlist(gup3)
  isnull <- gup3$short.desc=="NULL"
  strwrap.vec <- function(x){
    unlist(lapply(lapply(x, strwrap, width=strwrap.width), paste, collapse="\n"))
  }
  if(any(isnull)){
    gup3$short.desc[isnull] <- strwrap.vec(gup3$description[isnull])
  }
  gup3
}

removeOverlapsWithNullShortDesc <- function(gup, strwrap.width){
  genes <- unique(gup$hugo)
  up.list <- vector("list", length(genes))
  for(i in seq_along(genes)){
    up.list[[i]] <- reduceProtein(gup, genes[i], strwrap.width)
  }
  up.list <- GRangesList(up.list)
  up <- unlist(up.list)
  up
}


#' Extract features from uniprot database
#'
#' @param up uniprot features
#' @param fusions table of fusions with gene.5prime and gene.3prime columns
#' @param strwrap.width integer indicating how to wrap text for highly descriptive features in uniprot
#' @export
#' @examples
#' ## See fusions vignette
uniprotFeatures <- function(up, fusions, strwrap.width=30){
  up$feature.id <- seq_len(nrow(up))
  ## select features to plot
  feature_keys <- unique(up$feature_key)
  feature_keys <- feature_keys[-grep("residue", feature_keys)]
  feature_keys <- feature_keys[-grep("Repeat", feature_keys)]
  feature_keys <- feature_keys[c(1:3, 5, 7, 10, 11, 14, 16, 18, 20, 21, 22, 23, 24, 26, 27, 28, 30)]
  up2 <- up[up$feature_key %in% feature_keys, ]
  ## add back features that might have something to do with fusions
  up.fusion <- up[grep("fusion", up$description), ]
  up3 <- rbind(up2, up.fusion)
  up3 <- up3[!duplicated(up3$feature.id), ]
  up3 <- up3[order(up3$feature.id, decreasing=FALSE), ]

  ## some features overlap
  gup <- GRanges(up3$chrom, IRanges(start=up3$start, end=up3$stop),
                 description=up3$description,
                 shrt_desc=up3$shrt_desc,
                 hugo=up3$hugo,
                 aa_len=up3$aa_len)
  gup <- gup[gup$hugo %in% fusions$gene.5prime | gup$hugo %in% fusions$gene.3prime]
  gup2 <- removeOverlapsWithNullShortDesc(gup, strwrap.width)
  up2 <- as.data.frame(gup2)
  up2$end <- up2$end
  up2
}

bothGenesInUniprot <- function(x, uniprot){
  gene1.in <- x$gene.5prime %in% uniprot$hugo
  gene2.in <- x$gene.3prime %in% uniprot$hugo
  gene1.in & gene2.in
}

whichGeneInUniprot <- function(x, uniprot){
  gene1.in <- x$gene.5prime %in% uniprot$hugo
  gene2.in <- x$gene.3prime %in% uniprot$hugo
  setNames(c(gene1.in, gene2.in), c("gene1", "gene2"))
}

selectTx <- function(transcripts, fusions){
  genes <- c(fusions$gene.5prime, fusions$gene.3prime)
  ##tx.ids <- strsplit(fusions$fusion, "::")[[1]]
  tx.ids <- c(fusions$tx.5prime, fusions$tx.3prime)
  tx.ids <- gsub("\\(promoter\\)", "", tx.ids)
  roi <- transcripts[transcripts$tx_name %in% tx.ids]
  names(roi) <- roi$gene_name
  roi[genes]
}

genes <- function(x) {
  ##x <- unlist(strsplit(x$fusion, "::"))
  x <- c(x$tx.5prime, x$tx.3prime)
  gsub("\\(promoter\\)", "", x)
}

## copied from svovarian
exonsByTx <- function(ex, fusions){
  tx <- genes(fusions)
  tx1 <- tx[1]
  tx2 <- tx[2]
  ex1 <- ex[[tx1]]
  ex2 <- ex[[tx2]]
  tx1 <- GRanges(seqnames(ex1)[1],
                 IRanges(min(start(ex1)),
                         max(end(ex1))),
                 strand=strand(ex1)[1])
  tx2 <- GRanges(seqnames(ex2)[1],
                 IRanges(min(start(ex2)),
                         max(end(ex2))),
                 strand=strand(ex2)[1])
  tx <- c(tx1, tx2)
  tx$genes <- c(fusions$gene.5prime[1], fusions$gene.3prime[1])
  result <- list(exons1=ex1, exons2=ex2, transcripts=tx)
  names(result) <- c(tx$genes, "transcripts")
  result
}

meltReadPairs2 <- function(rdat, transcripts){
  gpairs <- improper(rdat)
  r1 <- as(first(gpairs), "GRanges")
  r2 <- as(GenomicAlignments::last(gpairs), "GRanges")
  names(r1) <- names(r2) <- NULL
  genes <- transcripts$genes
  r1$transcript <- genes[findOverlaps(r1, transcripts,
                                      select="first",
                                      ignore.strand=TRUE,
                                      maxgap=5000)]
  r2$transcript <- genes[findOverlaps(r2, transcripts,
                                      select="first",
                                      ignore.strand=TRUE,
                                      maxgap=5000)]
  r1$read <- "R1"
  r2$read <- "R2"
  r1$pair.id <- seq_len(length(r1))
  r2$pair.id <- r1$pair.id
  reads <- c(r1, r2)
  reads
}

.select_junction <- function(x){
  tab <- table(x)
  prop <- max(tab)/sum(tab)
  if(prop > 0.75){
    jxn <- as.integer(names(tab)[which.max(tab)])
  } else{
    if(diff(range(x)) < 5) {
      jxn <- median(jxn)
    } else {
      stop("no junction found")
    }
  }
  jxn
}

basepairJunction <- function(rreads, roi){
  split.reads <- rreads[rreads$is.split]
  gene1.reads <- split.reads[split.reads$transcript==names(roi)[1]]
  gene2.reads <- split.reads[split.reads$transcript==names(roi)[2]]
  stopifnot(identical(gene1.reads$pair.id, gene2.reads$pair.id))
  strands <- as.character(strand(roi))
  strand1 <- strands[1]
  strand2 <- strands[2]
  gene1.reads$strand_refgen <- strands[1]
  gene2.reads$strand_refgen <- strands[2]
  ## For tumor genome, assume coordinates of gene1 are same as in reference
  red.gene2 <- reduce(gene2.reads) 
  improper.reads <- rreads[!rreads$is.split]
  gene1.improper <- improper.reads[improper.reads$transcript==names(roi)[1]]
  gene2.improper <- improper.reads[improper.reads$transcript==names(roi)[2]]
  gene1.improper$strand_refgen <- strands[1]
  gene2.improper$strand_refgen <- strands[2]
  ## genome and gene2 has altered coordinates
  ##
  ## Always true: the end of the gene1 split reads (5' -> 3') abuts the
  ## start of the gene2 split reads 
  ## -- what changes is how we determine the start and end of the reads
  if(strand1=="-"){
    ##
    ## the end of the gene1 split reads (5' -> 3') abuts the start of
    ## the gene2 split reads 
    ##
    ## Since on minus strand, 'start' is the end of the gene1 split read
    bp.gene1 <- .select_junction(start(gene1.reads))
    ## trim improper reads that overhang the sequence junction
    starts <- pmax(bp.gene1, start(gene1.improper))
    start(gene1.improper) <- pmin(starts, end(gene1.improper))
  } else {
    bp.gene1 <- .select_junction(end(gene1.reads))
    ends <- pmin(bp.gene1[1], end(gene1.improper))
    end(gene1.improper) <- pmax(ends, start(gene1.improper))
  }
  if(strand2=="-"){
    ## Since on minus strand, 'end' is the start of the gene 2 split read
    bp.gene2 <- .select_junction(end(gene2.reads))
    ends <- pmin(bp.gene2[1], end(gene2.improper))
    end(gene2.improper) <- pmax(ends, start(gene2.improper))
  } else{
    bp.gene2 <- .select_junction(start(red.gene2))
    starts <- pmax(bp.gene2[1], start(gene2.improper))
    start(gene2.improper) <- pmin(starts, end(gene2.improper))
  }
  basepair.junction <- setNames(c(bp.gene1, bp.gene2), names(roi))
  rreads <- c(gene1.improper, gene2.improper, gene1.reads, gene2.reads)
  list(junction=basepair.junction,
       rearranged.reads=rreads)
}


# copied from svovarian
harmonizeReadMcols <- function(r.reads, reads){
  r.reads$pair.id <- r.reads$qname
  r.reads$pair.id <- as.integer(factor(r.reads$pair.id,
                                       levels=unique(r.reads$pair.id))) +
    max(reads$pair.id)
  ## combine the split reads with the aberrantly spaced reads
  reads2 <- granges(reads)
  ##reads2$rpid <- reads$rpid
  reads2$transcript <- reads$transcript
  reads2$read <- reads$read
  reads2$pair.id <- reads$pair.id
  reads2$is.split <- FALSE
  reads2$qname <- NA

  r.reads2 <- granges(r.reads)
  ##r.reads2$rpid <- r.reads$qname
  r.reads2$transcript <- r.reads$transcript
  r.reads2$read <- ""
  r.reads2$pair.id <- r.reads$pair.id
  r.reads2$is.split <- TRUE
  r.reads2$qname <- r.reads$qname
  reads3 <- c(reads2, r.reads2)
  reads3
}

exonTracks2 <- function(exons, reads, roi){
  genes <- names(roi)
  gene1 <- genes[1]
  gene2 <- genes[2]

  exons1 <- exons[[gene1]]
  exons2 <- exons[[gene2]]
  exons1$gene <- gene1
  exons2$gene <- gene2

  strands <- as.character(strand(roi))
  strand1 <- strands[1]
  strand2 <- strands[2]

  exons1$is_clipped <- NA
  exons2$is_clipped <- NA

  jxn1 <- as.integer(strsplit(roi$bp.jxn[1], ":")[[1]][[2]])
  jxn2 <- as.integer(strsplit(roi$bp.jxn[2], ":")[[1]][[2]])
  if(strand1=="-"){
    exons1$is_clipped <- ifelse(end(exons1) < jxn1, TRUE, FALSE)
  } else{
    exons1$is_clipped <- ifelse(start(exons1) > jxn1, TRUE, FALSE)
  }
  if(strand2=="-"){
    exons2$is_clipped <- ifelse(start(exons2) > jxn2, TRUE, FALSE)
  } else{
    exons2$is_clipped <- ifelse(end(exons2) < jxn2, TRUE, FALSE)
  }
  exons[[gene1]] <- exons1
  exons[[gene2]] <- exons2
  exons <- c(exons[[gene1]], exons[[gene2]])
  exons.df <- as(exons, "data.frame")
  ##mexons <- meltExons(exons)
  ##exons.df <- as(mexons, "data.frame")
  ##exons.df2 <- rbind(exons.df, clipped.df)
  e1 <- exons.df[exons.df$gene==genes[1], ]
  e2 <- exons.df[exons.df$gene==genes[2], ]
  exon_tracks <- list(gene1=e1, gene2=e2)
  names(exon_tracks) <- genes
  exon_tracks
}

.list_fusion_data <- function(fusion, exon_tracks, rearranged.reads){
  rreads <- rearranged.reads
  gene1 <- fusion$gene1
  gene2 <- fusion$gene2
  fusion_nm <- paste0(gene1, "-", gene2)
  id <- fusion$id
  id.rds <- paste0(id, ".rds")

  ##fusions <- readRDS(file.path("structuralvar/data/fusions/0fusions", id.rds))
  ##fusions <- fusions[grep(gene1, fusions$gene1), ]

  rreads <- reduce(rreads, ignore.strand=TRUE, min.gapwidth=100)
  names(rreads) <- unique(rearranged.reads$transcript)

  e1 <- exon_tracks[[1]]
  e2 <- exon_tracks[[2]]
  x2 <- max(e2$end)
  e2$midx <- apply(cbind(e2$start, e2$end), 1, mean)

  x <- min(e1$start)
  xend <- x+10e3
  e1$start <- e1$start-200
  e1$end <- e1$end+200
  e1$midx <- apply(cbind(e1$start, e1$end), 1, mean)
  exon_tracks[[1]] <- e1
  exon_tracks[[2]] <- e2
  e1$gene <- gene1
  e2$gene <- gene2
  names(exon_tracks)[1:2] <- c(gene1, gene2)
  list(rearranged.reads=rreads, exons=exon_tracks, fusion=fusion_nm)
}

fuseExonTracks <- function(data.list, basepair.jxn, roi){
  ##genes <- strsplit(data.list[["fusion"]], "-")[[1]]
  genes <- names(roi)
  gene1 <- genes[1]
  gene2 <- genes[2]
  tx1 <- data.list$exons[[gene1]]
  tx2 <- data.list$exons[[gene2]]
  tx1 <- tx1[!tx1$is_clipped, ]
  tx2 <- tx2[!tx2$is_clipped, ]
  strands <- as.character(strand(roi))
  ##
  ## assume gene1 is fixed (reference genome coordinates are correct)
  ##
  ## We must correct the coordinates for gene2
  basepair.jxn <- sapply(strsplit(basepair.jxn, ":"), "[", 2)
  basepair.jxn <- as.integer(basepair.jxn)
  starts <- tx2$start - basepair.jxn[2] + basepair.jxn[1]
  ends <- tx2$end - basepair.jxn[2] + basepair.jxn[1]
  tx2$start <- starts
  tx2$end <- ends
  ##
  ## tx1 is length-0 if promoter
  ##
  if(nrow(tx1) > 0){
    tx1$sequence_junction <- basepair.jxn[1]
  }
  tx2$sequence_junction <- basepair.jxn[2]
  #tx2$midx <- rowMeans(cbind(tx2$start, tx2$end))
  fused <- rbind(tx1, tx2)
  list(exons=fused, rearranged.reads=splitReads(data.list[["rlist"]]))
}

fuseExonTracksPM <- function(data.list, basepair.jxn, roi){
  ##genes <- strsplit(data.list[["fusion"]], "-")[[1]]
  genes <- names(roi)
  gene1 <- 1
  gene2 <- 2
  tx1 <- data.list$exons[[gene1]]
  tx2 <- data.list$exons[[gene2]]

  tx1 <- tx1[!tx1$is_clipped, ]
  tx2 <- tx2[!tx2$is_clipped, ]
  strands <- as.character(strand(roi))
  ##
  ## assume gene1 is fixed (reference genome coordinates are correct)
  ##
  ## We must correct the coordinates for gene2
  starts <- abs(tx2$start - basepair.jxn[gene2]) + basepair.jxn[gene1]
  ends <- abs(tx2$end - basepair.jxn[gene2]) + basepair.jxn[gene1]
  tx2$start <- starts
  tx2$end <- ends
  tx1$sequence_junction <- basepair.jxn[gene1]
  tx2$sequence_junction <- basepair.jxn[gene2]
  tx2$midx <- rowMeans(cbind(tx2$start, tx2$end))
  fused <- rbind(tx1, tx2)
  list(exons=fused, rearranged.reads=data.list$rearranged.reads)
}

proteinFeatures <- function(up, gene){
  features <- up[up$hugo==gene, ]
  features$midx <- rowMeans(cbind(features$start, features$end))
  features
}

proteinParams <- function(gene, is.first, description.size=1){
  if(is.first){
    p <- list(protein=gene,
              is.first=TRUE,
              background.color="lightblue",
              domain.color="steelblue",
              clipped.color="gray60")
  } else{
    p <- list(protein=gene,
              is.first=FALSE,
              background.color="beige",
              domain.color="darkolivegreen",
              clipped.color="gray60")
  }
  p$description.size <- description.size
  p
}

aaJunction <- function(rear, roi, cds.roi, bs.genome){
  tx.names <- setNames(roi$gene_name, roi$tx_name)
  cds.full <- fullTranscripts(cds.roi)
  names(cds.full) <- tx.names[names(cds.full)]
  cds.clipped <- clip(cds.roi)
  fused.tx <- fuse(cds.clipped)
  fused.protein <- tumorProtein(bs.genome, fused.tx)
  ref.protein <- referenceProtein(bs.genome, cds.full, tx.names)
  if(length(left(cds.clipped)) > 0){
    cds.clipped <- GRangesList(list(left(cds.clipped)[[1]], right(cds.clipped)[[1]]))
  } else{
    cds.clipped <- GRangesList(list(unlist(left(cds.clipped)),
                                    right(cds.clipped)[[1]]))
  }
  names(cds.clipped) <- tx.names[names(tx.names)]
  clipped.5prime <- referenceProtein(bs.genome, cds.clipped, tx.names[1])
  AA.break.5prime <- length(clipped.5prime[[1]])
  genes <- as.character(tx.names)
  nclipped.bases <- sum(width(cds.full[[genes[2]]])) -
    sum(width(cds.clipped[[genes[2]]]))
  nclipped.aa <- nclipped.bases/3
  AA.break.3prime <- nclipped.aa
  breaks <- c(AA.break.5prime, AA.break.3prime)
  breaks <- as.integer(round(breaks, 0))
  names(breaks) <- c("5'", "3'")
  breaks
}

clippedProtein <- function(p.dat, params, aa.jxn){
  if(params$is.first){
    clip.start <- aa.jxn
    clip.end <- p.dat$aa_len[1]
    p.dat$is.clipped <- p.dat$start > aa.jxn
  } else{
    clip.start <- 0
    clip.end <- aa.jxn
    p.dat$is.clipped <- p.dat$start < aa.jxn
  }
  p.dat$clip.start <- clip.start
  p.dat$clip.end <- clip.end
  p.dat$aa.jxn <- aa.jxn
  p.dat
}

fuseProteinTracks <- function(data.list){
  ## no strand issues with proteins -- already taken care of
  fusion <- rbind(data.list[[1]], data.list[[2]])
  ## add row that will the size of the clipped protein
  p1.domains <- data.list[[1]]
  ## coordinates of first protein 1 - aa.jxn
  p1 <- p1.domains[1, ]
  p1$start <- 1
  p1$end <- p1$aa.jxn
  ##p1$is.clipped <- FALSE
  ##p1$aa_len <- p1$end
  p1 <- data.frame(seqnames=p1$seqnames, start=p1$start, end=p1$end,
                   hugo=p1$hugo)
  aa.jxn <- p1$end
  junction.in.domain <- p1.domains$start < aa.jxn & aa.jxn < p1.domains$end
  if(any(junction.in.domain)){
    index <- which(junction.in.domain)
    p1.domains$end[index] <- p1.domains$end[index] <- aa.jxn
    p1.domains$midx <- rowMeans(cbind(p1.domains$start, p1.domains$end))
  }
  p1.domains <- p1.domains[p1.domains$end <= aa.jxn, ]

  ## p2  coordinates:  aa.jxn -> aa_len
  p2 <- data.list[[2]][1, ]
  aa.jxn2 <- p2$aa.jxn[1]
  p2$is.clipped <- FALSE
  p2$start <- aa.jxn
  p2$end <- p2$aa_len - aa.jxn2 + aa.jxn
  ##p2$aa_len <- p2$end
  p2 <- data.frame(seqnames=p2$seqnames, start=p2$start, end=p2$end,
                   hugo=p2$hugo)

  p2.domains <- data.list[[2]]
  ## does new junction occur within a domain
  junction.in.domain <- p2.domains$start < aa.jxn2 & aa.jxn2 < p2.domains$end
  if(any(junction.in.domain)){
    index <- which(junction.in.domain)
    p2.domains$end[index] <- p2.domains$end[index] <- aa.jxn2
  }
  p2.domains <- p2.domains[p2.domains$start >= aa.jxn2, , drop=FALSE]
  if(nrow(p2.domains) > 0){
    ## adjust start and end coordinates
    p2.domains$start <- p2.domains$start - aa.jxn2 + aa.jxn
    p2.domains$end <- p2.domains$end - aa.jxn2 + aa.jxn
    p2.domains$end <- p2.domains$end
    ##p2.domains$aa_len <- p2$aa_len
    p2.domains$aa.jxn <- aa.jxn
    p2.domains$midx <- rowMeans(cbind(p2.domains$start, p2.domains$end))
  }
  coords <- rbind(p1[1, ], p2[1, ])
  domains <- rbind(p1.domains, p2.domains)
  domains$aa.jxn <- rep(aa.jxn, nrow(domains))
  list(coords=coords, domains=domains)
}

geneParams <- function(roi){
  genes <- names(roi)
  strand <- as.character(strand(roi))
  gene1.params <- list(gene.name=genes[1],
                       background.color="lightblue",
                       exon.color="steelblue",
                       clipped.color="gray60",
                       is.first=TRUE,
                       strand=strand[1],
                       rearrangedStrand=strand[1])
  ##
  ## By convention, assume that the strand of gene 2 in the rearranged genome is
  ## on the same strand as gene 1.
  ##
  gene2.params <- list(gene.name=genes[2],
                       background.color="beige",
                       exon.color="darkolivegreen",
                       clipped.color="gray60",
                       is.first=FALSE,
                       strand=strand[2],
                       rearrangedStrand=strand[1])
  list(gene1.params, gene2.params)
}

## #' Collect transcript and protein-level data supporting a fusion for into a list object for subsequent plotting
## #'
## #'
## #' @param build length-one character vector providing genome build
## #' @param data.list a list object created by \code{listFusionData}
## #' @export
## #' @return a named list
## fusionData <- function(data.list, build=c("hg19", "hg18")){
##   x <- data.list[["fusions"]]
##   rlist <- data.list[["rlist"]]
##   uniprot <- data.list[["uniprot"]]
## 
##   transcripts <- loadTx(build)
##   roi <- selectTx(transcripts, x)
##   chroms <- as.character(seqnames(roi))
## 
##   rid <- x$rearrangement.id
##   if(length(rid) > 1) stop("no rearrangements")
##   r <- rlist[[rid]]
## 
##   txdb <- loadTxDb(match.arg(build))
##   cds.all <-   suppressWarnings(cdsBy(txdb, "tx",
##                                       use.names=TRUE))
##   exons.txdb <- suppressWarnings(exonsBy(txdb, "tx", use.names=TRUE))
## 
##   genes <- c(x$gene1, x$gene2)
##   exons <- exonsByTx(exons.txdb, x)
## 
##   reads <- meltReadPairs2(r, exons$transcripts)
##   readclusters <- exons$transcripts
## 
##   build <- "hg19"
##   bs.pkg <- paste0("BSgenome.Hsapiens.UCSC.", build)
##   bs.genome <- getBSgenome(bs.pkg)
## 
##   reads <- meltReadPairs2(r, exons$transcripts)
##   readclusters <- exons$transcripts
## 
##   r.reads <- sr[[rid]]
##   ix <- findOverlaps(r.reads, readclusters, select="first", maxgap=5000)
##   r.reads$transcript <- readclusters$genes[ix]
##   rreads <- harmonizeReadMcols(r.reads, reads)
##   bp.jxn <- basepairJunction(rreads, roi)
##   if(is.null(bp.jxn)){
##     reverse.roi <- TRUE
##     bp.jxn <- basepairJunction(rreads, roi[2:1])
##   } else reverse.roi <- FALSE
## 
##   rreads <- bp.jxn[["rearranged.reads"]]
##   bp.jxn <- bp.jxn[["junction"]]
##   roi$bp.jxn <- bp.jxn[names(roi)]
##   exon.tracks <- exonTracks2(exons, rreads, roi)
##   fusion.dat <- .list_fusion_data(x, exon.tracks, rreads)
## 
##   strands <- as.character(strand(roi))
##   strands <- paste(strands, collapse="")
##   fused.transcripts <- fuseExonTracks(fusion.dat, bp.jxn, roi)
##   if(strands=="+-"){
##     fused.transcripts <- fuseExonTracksPM(fusion.dat, bp.jxn, roi)
##   }
## 
##   p1 <- proteinFeatures(uniprot, genes[1])
##   p2 <- proteinFeatures(uniprot, genes[2])
##   p1.params <- proteinParams(genes[1], TRUE, description.size=2)
##   p2.params <- proteinParams(genes[2], FALSE, description.size=2)
##   g.params <- geneParams(roi)
## 
##   cds.roi <- getCDS(r, roi, cds.all)
##   aa.jxns <- aaJunction(r, roi, cds.roi, bs.genome)
##   roi$aa.jxn <- aa.jxns
## 
##   p1.clipped <- clippedProtein(p1, p1.params, aa.jxns["5'"])
##   p2.clipped <- clippedProtein(p2, p2.params, aa.jxns["3'"])
##   clip.list <- list(p1.clipped, p2.clipped)
##   names(clip.list) <- genes
##   p.fusions <- fuseProteinTracks(clip.list)
## 
##   rreads.df <- as.data.frame(rreads)
##   list(chroms=chroms,
##        roi=roi,
##        fusion.dat=fusion.dat,
##        g.params=g.params,
##        fused.transcripts=fused.transcripts,
##        rreads=rreads.df,
##        reverse.roi=reverse.roi,
##        protein1=p1,
##        protein2=p2,
##        protein1.clipped=p1.clipped,
##        protein2.clipped=p2.clipped,
##        protein1.params=p1.params,
##        protein2.params=p2.params,
##        protein.fusion=p.fusions)
## }


## fusionData2 <- function(data.list, build=c("hg19", "hg18")){
##   x <- data.list[["fusions"]]
##   rlist <- data.list[["rlist"]]
##   uniprot <- data.list[["uniprot"]]
## 
##   transcripts <- loadTx(build)
##   roi <- selectTx(transcripts, x)
##   chroms <- as.character(seqnames(roi))
## 
##   rid <- x$rearrangement.id
##   if(length(rid) > 1) stop("no rearrangements")
##   r <- rlist[[rid]]
## 
##   txdb <- loadTxDb(match.arg(build))
##   cds.all <-   suppressWarnings(cdsBy(txdb, "tx",
##                                       use.names=TRUE))
##   exons.txdb <- suppressWarnings(exonsBy(txdb, "tx", use.names=TRUE))
## 
##   genes <- c(x$gene1, x$gene2)
##   exons <- exonsByTx(exons.txdb, x)
## 
##   reads <- meltReadPairs2(r, exons$transcripts)
##   readclusters <- exons$transcripts
## 
##   build <- "hg19"
##   bs.pkg <- paste0("BSgenome.Hsapiens.UCSC.", build)
##   bs.genome <- getBSgenome(bs.pkg)
## 
##   reads <- meltReadPairs2(r, exons$transcripts)
##   readclusters <- exons$transcripts
## 
##   r.reads <- sr[[rid]]
##   ix <- findOverlaps(r.reads, readclusters, select="first", maxgap=5000)
##   r.reads$transcript <- readclusters$genes[ix]
##   rreads <- harmonizeReadMcols(r.reads, reads)
##   bp.jxn <- basepairJunction(rreads, roi)
##   if(is.null(bp.jxn)){
##     reverse.roi <- TRUE
##     bp.jxn <- basepairJunction(rreads, roi[2:1])
##   } else reverse.roi <- FALSE
## 
##   rreads <- bp.jxn[["rearranged.reads"]]
##   bp.jxn <- bp.jxn[["junction"]]
##   roi$bp.jxn <- bp.jxn[names(roi)]
##   exon.tracks <- exonTracks2(exons, rreads, roi)
##   fusion.dat <- .list_fusion_data(x, exon.tracks, rreads)
## 
##   strands <- as.character(strand(roi))
##   strands <- paste(strands, collapse="")
##   fused.transcripts <- fuseExonTracks(fusion.dat, bp.jxn, roi)
##   if(strands=="+-"){
##     fused.transcripts <- fuseExonTracksPM(fusion.dat, bp.jxn, roi)
##   }
## 
##   p1 <- proteinFeatures(uniprot, genes[1])
##   p2 <- proteinFeatures(uniprot, genes[2])
##   p1.params <- proteinParams(genes[1], TRUE, description.size=2)
##   p2.params <- proteinParams(genes[2], FALSE, description.size=2)
##   g.params <- geneParams(roi)
## 
##   cds.roi <- getCDS(r, roi, cds.all)
##   aa.jxns <- aaJunction(r, roi, cds.roi, bs.genome)
##   roi$aa.jxn <- aa.jxns
## 
##   p1.clipped <- clippedProtein(p1, p1.params, aa.jxns["5'"])
##   p2.clipped <- clippedProtein(p2, p2.params, aa.jxns["3'"])
##   clip.list <- list(p1.clipped, p2.clipped)
##   names(clip.list) <- genes
##   p.fusions <- fuseProteinTracks(clip.list)
## 
##   rreads.df <- as.data.frame(rreads)
##   list(chroms=chroms,
##        roi=roi,
##        fusion.dat=fusion.dat,
##        g.params=g.params,
##        fused.transcripts=fused.transcripts,
##        rreads=rreads.df,
##        reverse.roi=reverse.roi,
##        protein1=p1,
##        protein2=p2,
##        protein1.clipped=p1.clipped,
##        protein2.clipped=p2.clipped,
##        protein1.params=p1.params,
##        protein2.params=p2.params,
##        protein.fusion=p.fusions)
## }

#' Collect fusion-related data into a single list
#'
#' @param fusions an object returned by \code{fusionList}
#' @param rlist a \code{RearrangementList}
#' @return a named list 
#' @export
listFusionData <- function(rlist, fusions){
  extdata <- system.file("extdata", package="svfusions")
  up <- readRDS(file.path(extdata, "uniprot.rds"))
  up2 <- uniprotFeatures(up, fusions, strwrap.width=20)
  both.genes <- bothGenesInUniprot(fusions, up2)
  if(!any(both.genes)){
    in.uniprot <- whichGeneInUniprot(fusions, up2)
    msg <- "Both genes are not in uniprot."
    message(msg)
    return(in.uniprot)
  }
  x <- fusions[both.genes, ]
  List(rlist=rlist,
       fusions=x,
       uniprot=up2,
       split_reads=splitReads(rlist))
}
