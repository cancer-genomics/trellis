context("noncodingRanges")

test_that("noncoding placeholders match either transcripts schema", {
  g <- GRanges(c("chr1", "chr2"), IRanges(c(100, 200), width = 10),
               extra = c("a", "b"))
  old <- GRanges("chr1", IRanges(1, 10), tx_id = 1L, tx_name = "NM_1",
                 gene_name = "A", cancer_connection = TRUE, biol_sign = TRUE)
  both <- old
  both$clinically_significant <- "Tx1"
  both$cancer_gene <- TRUE
  for (tx in list(old, both)) {
    nc <- trellis:::noncodingRanges(g, tx)
    expect_identical(names(mcols(nc)), names(mcols(tx)))
    expect_identical(vapply(as.list(mcols(nc)), class, ""),
                     vapply(as.list(mcols(tx)), class, ""))
    expect_identical(nc$gene_name, c("noncoding1", "noncoding2"))
    expect_identical(nc$tx_name, c("", ""))
    expect_length(c(tx, nc), 3L)
  }
  nc <- trellis:::noncodingRanges(g, both)
  ## NA, not FALSE: driver_genes(clin_sign = TRUE) selects on !is.na().
  expect_true(all(is.na(nc$clinically_significant)))
  expect_false(any(nc$cancer_gene))
})
