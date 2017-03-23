context("Revising sequence junctions")

.test_that <- function(name, expression) NULL

## overlapping hemizygous deletions on chromosome 3

cgov10t_preprocess <- function(){
  chroms <- c("chr3", "chr3", "chr3", "chr4", "chr9", "chrX")
  starts <- c(59600000, 60246000, 175025000, 12660000, 9800000, 6686000)
  ends <- c(61000000, 60700000, 175167000, 12770000, 10500000, 7000000)
  bed <- GRanges(chroms, IRanges(starts, ends))
  bed <- keepSeqlevels(bed, "chr3", pruning.mode="coarse")

  bamdir <- system.file("extdata", package="svbams", mustWork=TRUE)
  bamfile <- file.path(bamdir, "cgov10t.bam")
  svfile <- file.path(bamdir, "cgov10t_deletions.rds")
  expected.sv <- readRDS(svfile)
  expected.irp <- improper(expected.sv)

  cbsfile <- file.path(bamdir, "cgov10t_cbs.rds")
  lrfile <- file.path(bamdir, "cgov10t_binlevel.rds")
  bins <- readRDS(file.path(bamdir, "cgov10t_bins1kb.rds"))
  ##
  ## todo: regenerate irp from bamfile
  ##
  ## make sure all improper read pairs are available
  ##irpfile <- file.path(bamdir, "cgov10t_irp.rds")
  ##irp <- readRDS(irpfile)
  lr <- readRDS(lrfile)
  bins$log_ratio <- lr/1000
  segments <- readRDS(cbsfile)
  segments <- subsetByOverlaps(segments, bed)
  read_pairs <- listReadPairs(bamfile, segments)
  pdat <- preprocessData(bam.file=bamfile,
                         genome=genome(bins)[[1]],
                         bins=bins,
                         segments=segments,
                         read_pairs=read_pairs)
}

test_that("overlappingHemizgyous", {
  library(Rsamtools)
  library(svbams)
  pdat <- cgov10t_preprocess()
  sv1 <- sv_deletions(pdat)
  if(FALSE){
    df <- data.frame(start=start(bins),
                     lr=bins$log_ratios/1000)
    library(ggplot2)
    ggplot(df, aes(start, lr)) + geom_point(size=0.5) +
      coord_cartesian(xlim=c(60244000, 60697000))
    irp1 <- improper(expected.sv[1])
    fst <- as(granges(first(irp1)), "data.frame")
    lst <- as(granges(last(irp1)), "data.frame")
    fst$start2 <- lst$start
    fst$end2 <- lst$end
    starts <- apply(fst[, c("start", "start2")], 1, min)
    ends <- apply(fst[, c("end", "end2")], 1, max)
    df2 <- data.frame(start=starts, end=ends, y=seq_along(starts))
    ggplot(df, aes(start, lr)) + geom_point(size=0.5) +
      coord_cartesian(xlim=c(60244000, 60697000))   +
      geom_segment(data=df2, aes(x=start, xend=end, y=y, yend=y),
                   inherit.aes=FALSE)
  }
  sv <- deletion_call(pdat)
  param <- DeletionParam()
  calls(sv) <- rpSupportedDeletions(sv, param=param, bins=preprocess$bins)
  sv <- removeHemizygous(sv)
  sv <- reviseEachJunction(sv, pdat$bins, pdat$improper_rp, param)
  ##pdat$improper_rp, param)
  sv.revise <- revise(sv, pdat$bins, param)
  if(FALSE){
    df3 <- data.frame(start=start(variant(sv)),
                      end=end(variant(sv)),
                      y=variant(sv)$seg.mean)
    ggplot(df, aes(start, lr)) + geom_point(size=0.5) +
      coord_cartesian(xlim=c(60244000, 60697000))   +
      geom_segment(data=df2, aes(x=start, xend=end, y=y, yend=y),
                   inherit.aes=FALSE) +
      geom_segment(data=df3, aes(x=start, xend=end, y=y, yend=y),
                   inherit.aes=FALSE)
  }
  copynumber(sv) <- granges_copynumber2(variant(sv), pdat$bins)
  calls(sv) <- rpSupportedDeletions(sv, param=param, pdat$bins)
  indexImproper(sv) <- updateImproperIndex(sv, maxgap=500)
  calls(sv) <- rpSupportedDeletions(sv, param, pdat$bins)
  sv2 <- leftHemizygousHomolog(sv, pdat$bins, param)
  sv3 <- rightHemizygousHomolog(sv2, pdat$bins, param)
  calls(sv3) <- rpSupportedDeletions(sv3, param, pdat$bins)
  message("Refining homozygous boundaries by spanning hemizygous+")
  sv5 <- refineHomozygousBoundaryByHemizygousPlus(sv3)
  sv6 <- callOverlappingHemizygous(sv5)
  sv7 <- removeSameStateOverlapping2(sv6)
  sv <- removeHemizygous(sv7)
  expect_identical(sv, sv.revise)

  sv2 <- rename(sort(sv))
  ##expected <- rename(sort(expected.sv))
  g <- granges(variant(sv2))
  ##g.expected <- granges(variant(expected))
  ##expect_identical(start(g), start(g.expected))
  ##expect_identical(end(g), end(g.expected))
  ##expect_identical(calls(sv2), calls(expected))
  gr_filters <- genomeFilters("hg19")
  sv.finalize <- finalize_deletions(sv, pdat, gr_filters, param)
  expect_identical(variant(sv.finalize),
                   variant(sv1))
  if(FALSE){
    df3 <- data.frame(start=start(variant(sv)),
                      end=end(variant(sv)),
                      y=variant(sv)$seg.mean)
    ggplot(df, aes(start, lr)) + geom_point(size=0.5) +
      coord_cartesian(xlim=c(60244000, 60697000))   +
      geom_segment(data=df2, aes(x=start, xend=end, y=y, yend=y),
                   inherit.aes=FALSE) +
      geom_segment(data=df3, aes(x=start, xend=end, y=y, yend=y),
                   inherit.aes=FALSE)

    df4 <- data.frame(start=start(variant(expected.sv)),
                      end=end(variant(expected.sv)),
                      y=copynumber(expected.sv))
    ggplot(df, aes(start, lr)) + geom_point(size=0.5) +
      coord_cartesian(xlim=c(60244000, 60697000))   +
      geom_segment(data=df2, aes(x=start, xend=end, y=y, yend=y),
                   inherit.aes=FALSE) +
      geom_segment(data=df4, aes(x=start, xend=end, y=y, yend=y),
                   inherit.aes=FALSE, color="blue")
  }
})

.test_that <- function(name, expr) NULL
.test_that("scratch", {
  g.old <- variant(sv.old)
  g.cur <- variant(sv)
  sv.delta <- sv[!overlapsAny(g.cur, g.old)]


  library(ggplot2)
  suppressPackageStartupMessages(library(gridExtra))
  library(scales)

  segs <- readRDS(file.path("structuralvar/data/segment/0cbs", id.rds))
  segs <- segs[segs$seg.mean < log2(0.75)]
  suppressPackageStartupMessages(library(SummarizedExperiment))
  ##
  ## region of interest (roi)
  ##
  roi <- variant(sv.delta[1])
  roi <- keepSeqlevels(roi, as.character(seqnames(roi)))
  expand <- 200e3
  roi2 <- GRanges(seqnames(roi), IRanges(start(roi)-expand,
                                         end(roi) + expand))
  seqinfo(roi2) <- seqinfo(roi)
  pview2 <- pview[queryHits(findOverlaps(rowRanges(pview), roi2)), ]
  segs <- keepSeqlevels(segs, seqlevels(roi), pruning.mode="coarse")
  segs.df <- as(segs, "data.frame")
  df <- data.frame(logr=assays(pview2)[, 1],
                   start=start(rowRanges(pview2)))
ylim <- c(-9, 2)
df$logr <- svpreprocess::threshold(df$logr, ylim)
brks <- pretty(df$start, n=8)
region <- subsetByOverlaps(segs, roi)
region <- region[region$seg.mean < -1]
region <- as.data.frame(region)
xlim <- c(start(roi2), end(roi2))

A <- ggplot(df, aes(start, logr)) +
  geom_point(size=1, color="gray50") +
  scale_x_continuous(expand=c(0,0), breaks=brks, labels=brks/1e6)+
  scale_y_continuous(expand=c(0,0)) +
  geom_segment(data=segs.df,
               aes(x=start, xend=end, y=seg.mean, yend=seg.mean),
               size=1) +
  coord_cartesian(xlim=xlim, ylim=ylim) +
  ylab(expression(log[2]~ratio)) +
  geom_rect(data=region,
            aes(xmin=start, xmax=end, ymin=-Inf, ymax=+Inf),
            fill="steelblue", color="transparent", alpha=0.3,
            inherit.aes=FALSE) +
  theme(axis.text=element_text(size=10),
        axis.text.x=element_blank()) + xlab("") +
  annotate("text", x=xlim[1] + 15e3, y=-8, label="chr15", size=3)
A1 <- ggplotGrob(A)

  rps <- thinReadPairs(sv.delta[1])
  library(svfilters.hg19)
  ##trace(svcnvs:::meltReadPairs, browser)
  rps <- svcnvs:::meltReadPairs(rps)
colors <- c("#0072B2", "#009E73")
p <- ggplot(rps, aes(ymin=readpair-0.2, ymax=readpair+0.2,
                xmin=start/1e6, xmax=end/1e6, color=read,
                fill=read, group=readpair)) +
  geom_rect() +
  xlim(c(min(rps$start), max(rps$end))/1e6) +
  geom_line(aes(x=start/1e6, y=readpair)) +
  ylab("read pair index") +
  scale_x_continuous(breaks=pretty_breaks(5)) +
  geom_rect(data=region,
            aes(xmin=start/1e6, xmax=end/1e6, ymin=-Inf, ymax=+Inf),
            fill="steelblue", color="transparent", alpha=0.2,
            inherit.aes=FALSE) +
  scale_color_manual(values=colors) +
  scale_fill_manual(values=colors) + 
  xlab("Mb") +
  theme(axis.text.x=element_text(size=7)) +
  guides(fill=FALSE, color=FALSE) +
  geom_vline(xintercept=c(start(roi)/1e6, end(roi)/1e6), linetype="dashed")
  B <- ggplotGrob(p)
})
