context("driver_genes() schema tolerance (card EXT-01)")

## driver_genes() must resolve either of two annotation schemas on its
## transcripts GRanges:
##   - current schema:    'cancer_gene' / 'clinically_significant'
##   - historical schema: 'biol_sign' / 'cancer_connection' (what every
##     released svfilters.hg18/hg19 'transcripts' object actually ships)
## and must fail loudly -- never silently return an empty result -- when
## neither schema is present.

.mock_tx_current <- function(){
  GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = c(1, 1001, 2001, 3001), width = 100),
    gene_name = c("GENE_A", "GENE_B", "GENE_C", "GENE_D"),
    cancer_gene = c(TRUE, FALSE, TRUE, FALSE),
    clinically_significant = c("Level 1", NA, NA, NA)
  )
}

.mock_tx_historical <- function(){
  GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = c(1, 1001, 2001, 3001), width = 100),
    gene_name = c("GENE_A", "GENE_B", "GENE_C", "GENE_D"),
    biol_sign = c(TRUE, FALSE, TRUE, FALSE),
    cancer_connection = c(TRUE, FALSE, FALSE, FALSE)
  )
}

.mock_tx_unrecognized <- function(){
  GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = c(1, 1001), width = 100),
    gene_name = c("GENE_A", "GENE_B"),
    some_other_column = c(TRUE, FALSE)
  )
}

test_that("driver_genes() resolves the current schema", {
  tx <- .mock_tx_current()
  expect_identical(driver_genes(tx, clin_sign = FALSE), c("GENE_A", "GENE_C"))
  expect_identical(driver_genes(tx, clin_sign = TRUE), "GENE_A")
})

test_that("driver_genes() resolves the historical schema", {
  tx <- .mock_tx_historical()
  expect_identical(driver_genes(tx, clin_sign = FALSE), c("GENE_A", "GENE_C"))
  expect_identical(driver_genes(tx, clin_sign = TRUE), "GENE_A")
})

test_that("driver_genes() fails loudly -- not silently -- on an unrecognized schema", {
  tx <- .mock_tx_unrecognized()
  expect_error(driver_genes(tx, clin_sign = FALSE),
               "could not resolve a recognized schema")
  expect_error(driver_genes(tx, clin_sign = TRUE),
               "could not resolve a recognized schema")
  ## the error must report both what was found and what was expected, so a
  ## future reader isn't left guessing why the match failed
  err <- tryCatch(driver_genes(tx, clin_sign = FALSE), error = function(e) e)
  expect_match(conditionMessage(err), "some_other_column", fixed = TRUE)
  expect_match(conditionMessage(err), "cancer_gene", fixed = TRUE)
  expect_match(conditionMessage(err), "biol_sign", fixed = TRUE)
})

## recurrentDrivers() (amplicon-utils.R:1362) had the same class of defect --
## it read a column name ('cancer_genes', plural) that no writer in the
## package ever produces. setDrivers()/getDrivers()/standardizeGRangesMetadata()
## all write 'cancer_gene' (singular); the pre-2016 code wrote 'driver'.

.mock_amplicon_current <- function(){
  GRangesList(
    sample1 = GRanges(
      seqnames = "chr1",
      ranges = IRanges(start = c(1, 1001), width = 100),
      cancer_gene = c("GENE_A", NA)
    )
  )
}

.mock_amplicon_historical <- function(){
  GRangesList(
    sample1 = GRanges(
      seqnames = "chr1",
      ranges = IRanges(start = c(1, 1001), width = 100),
      driver = c("GENE_A", NA)
    )
  )
}

.mock_amplicon_unrecognized <- function(){
  GRangesList(
    sample1 = GRanges(
      seqnames = "chr1",
      ranges = IRanges(start = c(1, 1001), width = 100),
      some_other_column = c("GENE_A", NA)
    )
  )
}

.mock_gene_transcripts <- function(){
  GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = 1, width = 100),
    gene_name = "GENE_A"
  )
}

test_that("recurrentDrivers() resolves the current schema", {
  grl <- .mock_amplicon_current()
  out <- recurrentDrivers(grl, transcripts = .mock_gene_transcripts())
  expect_identical(out$gene, "GENE_A")
})

test_that("recurrentDrivers() resolves the historical schema", {
  grl <- .mock_amplicon_historical()
  out <- recurrentDrivers(grl, transcripts = .mock_gene_transcripts())
  expect_identical(out$gene, "GENE_A")
})

test_that("recurrentDrivers() fails loudly -- not silently -- on an unrecognized schema", {
  grl <- .mock_amplicon_unrecognized()
  expect_error(recurrentDrivers(grl, transcripts = .mock_gene_transcripts()),
               "could not resolve a recognized driver-gene column")
  err <- tryCatch(recurrentDrivers(grl, transcripts = .mock_gene_transcripts()),
                   error = function(e) e)
  expect_match(conditionMessage(err), "some_other_column", fixed = TRUE)
  expect_match(conditionMessage(err), "cancer_gene", fixed = TRUE)
  expect_match(conditionMessage(err), "driver", fixed = TRUE)
})
