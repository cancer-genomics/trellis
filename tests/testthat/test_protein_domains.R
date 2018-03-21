context("Protein domains")

test_that("protein_domains", {
  extdata <- system.file("extdata", package="svfusions")
  fusions <- readRDS(file.path(extdata, "valid_fusions.rds"))
  tab <- fusionTable2(fusions)
  up <- readRDS(file.path(extdata, "uniprot.rds"))
  up2 <- uniprotFeatures(up, tab, strwrap.width=20)
  ## check that none of the short descriptions are too long
  nchar.desc <- nchar(up2$short.desc)
  gene1.in.up <- tab$gene.5prime %in% up2$hugo
  gene2.in.up <- tab$gene.3prime %in% up2$hugo
  bothin.up <- gene1.in.up & gene2.in.up
  tab <- tab[bothin.up, ]
  ## the aa_len field of the amino acid in uniprot might have a zero-based index
  ## -- it is one smaller than the length we've recorded
  ## coerce to GRanges
  tumor_aa_ranges <- aa_granges(tab)
  expect_identical(as.character(seqnames(tumor_aa_ranges)),
                   rep(c("ERBB4", "IKZF2"), each=2))
  domain_aa_ranges <- GRanges(up2$hugo, IRanges(up2$start, up2$end),
                              chromosome=up2$seqnames,
                              description=up2$description,
                              short.desc=up2$short.desc,
                              aa_len=up2$aa_len)
  domains <- subsetByOverlaps(domain_aa_ranges, tumor_aa_ranges, type="within")
  unique(domains$short.desc)
  expect_identical(length(domains), 4L)
})

.test_that <- function(nm, expr) NULL

.test_that("", {
  ##
  ## Plots
  ##
  ## how to test the plotting functions?
  load_all("../../../svplots")
  g.ideo <- ggIdeogram(dat[["chroms"]],
                       dat[["fusion.dat"]],
                       dat[["g.params"]])
  g1.clip <- ggClippedExons(dat[["fusion.dat"]],
                            dat[["g.params"]][[1]],
                            dat[["roi"]])
  g.axis <- ggAxisLabel(dat[["fusion.dat"]],
                        dat[["g.params"]][[1]], "3'")
  gg1.clip <- ggplotGrob(g1.clip)
  gg.axis <- ggplotGrob(g.axis)
  gg.axis$heights <- gg1.clip$heights
  g2.clip <- ggClippedExons(dat[["fusion.dat"]],
                            dat[["g.params"]][[2]],
                            dat[["roi"]])
  gg2.clip <- ggplotGrob(g2.clip)
  g.fused <- ggTxFusion(dat[["fused.transcripts"]],
                        dat[["g.params"]][[1]],
                        dat[["g.params"]][[2]])
  g.fused <- g.fused + ggtitle(dat[["fusion.dat"]]$fusion)
  gg.fused <- ggplotGrob(g.fused)
  gg.readlist <- ggRearrangedReads(dat[["rreads"]], dat[["g.params"]],
                                   dat[["roi"]], dat[["reverse.roi"]])
  ## plot protein
  g.p1 <- ggProtein(dat[["protein1"]],
                    dat[["protein1.params"]])
  g.p2 <- ggProtein(dat[["protein2"]],
                    dat[["protein2.params"]])
  ##load_all("integration/integration")
  gp1.clip <- ggClippedProtein1(dat[["protein1.clipped"]],
                                dat[["protein1.params"]])
  gp2.clip <- ggClippedProtein2(dat[["protein2.clipped"]],
                                dat[["protein2.params"]])
  ##load_all("integration/integration")
  gp.fuse <- ggProteinFusion(dat[["protein.fusion"]],
                             dat[["protein1.params"]],
                             dat[["protein2.params"]])

  cp <- compositeFusionParams()
  blank <- gg_blank()
  id2 <- gsub("\\.bam", "", id)
  library(gridExtra)
  grid.arrange(g.ideo,
               gg1.clip, gg.axis,
               gg2.clip, gg.axis,
               gg.readlist[[1]],
               gg.readlist[[2]],
               gg.fused, gg.axis,
               g.p1, g.p2,
               gp1.clip, gp2.clip,
               blank, gp.fuse,
               widths=cp[["widths"]],
               heights=cp[["heights"]],
               layout_matrix=cp[["layout"]])
})



.test_that("protein domains", {
  data(rear_cgov7t, package="svfusions")
  extdata <- system.file("extdata", package="svfusions")
  x <- readRDS(file.path(extdata, "fusionList.641eb81.rds"))
  flist <- listFusionData(rlist=rear_cgov7t, fusions=x[1, ])
  dat <- fusionData(data.list=flist)
  if(FALSE){
    saveRDS(dat, file="fusionData.664a371.rds")
  }
  expected <- readRDS("fusionData.664a371.rds")
  expect_equivalent(expected, dat)
})
