context("ggplot for Rearrangement class")

test_that("ggRearrange", {
  extdata <- system.file("extdata", package="svrearrange")
  rfile <- file.path(extdata, "CGOV11T_1.bam.rds")
  rlist <- readRDS(rfile)
  r <- rlist[[1]]
  ##trace(rearDataFrame, browser)
  df <- rearDataFrame(r, "hg19")
  ##trace(ggRearrange, browser)
  p <- ggRearrange(df)
  expect_is(p, "gg")

  ## a more challenging example
  id <- "10546-10582"
  r <- rlist[[id]]
  df <- rearDataFrame(r, "hg19")
  if(FALSE)
    ggRearrange(df)
  expect_identical(levels(df$region),
                   c("C3orf67-AS1,C3orf67", "FHIT"))

  index <- c(1, match(id, names(rlist)))
  r.ordered <- fiveTo3List(rlist[index], "hg19")

  df <- rearDataFrameList(r.ordered[1:2], "hg19")
  if(FALSE) ggRearrange(df)
  df <- rearDataFrameList(r.ordered[3:4], "hg19")
  if(FALSE) ggRearrange(df)
})
