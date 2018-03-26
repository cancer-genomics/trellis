context("Protein domains")

test_that("protein_domains", {
  extdata <- system.file("extdata", package="trellis")
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


