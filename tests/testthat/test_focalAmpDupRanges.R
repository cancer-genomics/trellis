context("Duplication events surrounding focal amplicons")

test_that("focalAmpliconDupRanges", {
  params <- ampliconParams()
  path <- system.file("extdata", package="svbams")
  ag <- readRDS(file.path(path, "addFocalDups.ffab104.rds"))

  LOW_THR <- params[["LOW_THR"]]
  MAX_SIZE <- params[["maxgap"]]
  g <- ampliconRanges(ag)
  is_dup <- isDuplication(ranges(ag), minimum_foldchange=LOW_THR)
  tab <- as.numeric(table(is_dup))
  expect_identical(tab, c(224, 6))
  ## remove any lower copy ranges that are very large
  dup_g <- ranges(ag)[is_dup]
  dup_g <- dup_g[width(dup_g) < MAX_SIZE]
  g <- sort(c(dup_g, g))
  ##
  ## make the set of amplicon ranges and low-copy duplications bigger by 1kb (why?)
  ##
  if(length(g) > 0){
    g <- expandGRanges(g, 1e3)
  }
  ## we do not want to link germline events
  germ <- germline(ag)
  ##
  ## TODO: add expansion to params
  ##
  germ <- reduce(expandGRanges(germ, 10e3))
  g2 <- filterBy(g, germ, type="within")
  g3 <- g[!overlapsAny(g, germ, type="within")]
  expect_identical(g2, g3)
  expect_identical(length(g3), 10L)
  if(FALSE){
    segs <- as.data.frame(ranges(ag))
    segs <- segs[!is.na(segs$seg.mean), ]
    lowcopy.dups <- as(dup_g, "data.frame")
    ggplot(segs, aes(x=start,
                     xend=end,
                     y=seg.mean,
                     yend=seg.mean)) +
      geom_segment(aes(color=is_amplicon), size=3) +
      geom_segment(data=lowcopy.dups, color="orange", size=3) +
      scale_color_manual(values=c("gray", "blue")) +
      geom_hline(yintercept=params[["LOW_THR"]], linetype="dashed") +
      facet_wrap(~seqnames) +
      ylim(c(-3, 3)) +
      xlab("Mb")

    segs.chr8 <- segs[segs$seqnames=="chr8", ]
    lc.chr8 <- lowcopy.dups[lowcopy.dups$seqnames=="chr8", ]
    ggplot(segs.chr8, aes(x=start/1e6,
                     xend=end/1e6,
                     y=seg.mean,
                     yend=seg.mean)) +
      geom_segment(aes(color=is_amplicon), size=3) +
      geom_segment(data=lc.chr8, color="orange", size=3) +
      scale_color_manual(values=c("gray", "blue")) +
      geom_hline(yintercept=params[["LOW_THR"]], linetype="dashed") +
      ylim(c(-3, 3)) +
      coord_cartesian(xlim=c(127, 131)) +
      xlab("Mb")

    ggplot(segs.chr8, aes(x=start/1e6,
                     xend=end/1e6,
                     y=seg.mean,
                     yend=seg.mean)) +
      geom_segment(aes(color=is_amplicon), size=3) +
      geom_segment(data=lc.chr8, color="orange", size=3) +
      scale_color_manual(values=c("gray", "blue")) +
      geom_hline(yintercept=params[["LOW_THR"]], linetype="dashed") +
      ylim(c(-3, 3)) +
      coord_cartesian(xlim=c(140, 146)) +
      xlab("Mb")
  }
})
