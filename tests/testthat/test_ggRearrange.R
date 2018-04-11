context("ggplot for Rearrangement class")

test_that("ggRearrange", {
  extdata <- system.file("extdata", package="trellis")
  rfile <- file.path(extdata, "CGOV11T_1.bam.rds")
  rlist <- readRDS(rfile)
  r <- rlist[[1]]
  r2 <- fiveTo3Prime(r, "hg19")
  df <- rearDataFrame(r2[[1]], "hg19")
  p <- ggRearrange(df)
  expect_is(p, "list")

  ## a more challenging example
  id <- "10546-10582"
  r <- rlist[[id]]
  r <- fiveTo3Prime(r, "hg19")[[1]]
  df <- rearDataFrame(r, "hg19")
  if(FALSE)
    ggRearrange(df)
  expect_identical(levels(df$region),
                   c("C3orf67-AS1,C3orf67", "FHIT"))

  index <- match(id, names(rlist))
  r.ordered <- fiveTo3List(rlist[index], "hg19")

  df <- rearDataFrameList(r.ordered[1:2], "hg19")
  if(FALSE) ggRearrange(df)
})
