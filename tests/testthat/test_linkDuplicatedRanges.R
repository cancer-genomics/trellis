test_that("linkDuplicatedRanges", {
  left_first <- GRanges("chrY", IRanges(58967001, 59000001), seg.mean=1.07, overlaps_germline="no_germline")
  left_last <- GRanges("chrY", IRanges(59001001, 59032001), seg.mean=2.82, overlaps_germline="no_germline")
  right_first <- GRanges(rep("chrY", 2),
                         IRanges(c(13401001, 59001001),
                                 c(13412001, 59032001)),
                         seg.mean=c(2.29, 2.82), overlaps_germline=c("no_germline", "no_germline"))
  right_last <- GRanges("chrY", IRanges(13414001, 13442001), seg.mean=0.88, overlaps_germline="no_germline")
  left <- list(first=left_first, last=left_last)
  right <- list(first=right_first, last=right_last)
  flank <- list(left=left, right=right)
  gapsLeft <- GRanges(seqnames(flank[["left"]]$first),
                      IRanges(end(flank[["left"]]$first)+1,
                              start(flank[["left"]]$last)))
  ## An error occurs here because flank[["right"]]$first is a length-2 GRanges and
  ## flank[["right"]]$last is a length-1 GRanges.
  expect_error(gapsRight <- GRanges(seqnames(flank[["right"]]$first),
                                    IRanges(end(flank[["right"]]$first)+1,
                                            start(flank[["right"]]$last))))
  
})
