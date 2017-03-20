context("Revising sequence junctions")

.test_that <- function(name, expression) NULL

## overlapping hemizygous deletions on chromosome 3

.test_that("scratch", {
  library(Rsamtools)
  library(svbams)
  bamdir <- system.file("extdata", package="svbams")
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
  irpfile <- file.path(bamdir, "cgov10t_irp.rds")
  irp <- readRDS(irpfile)
  mapped <- irp[!is.na(as.character(seqnames(irp)))]
  expect_true(all(overlapsAny(expected.irp, mapped, type="within")))

  gr <- readRDS(cbsfile)
  gr <- gr[gr$seg.mean < log2(0.75)]

  bview <- BamViews(bamPaths=bamfile, bamRanges=bins)
  aview <- AlignmentViews2(bview, irpfile)
  pview <- PreprocessViews2(bview)
  paths(pview) <- lrfile
  setScale(pview) <- 1000

  set.seed(123)
  sv1 <- sv_deletions(gr=gr,
                      aview=aview,
                      bview=bview,
                      pview=pview)
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

  sv <- deletion_call(aview, pview, gr)
  calls(sv) <- rpSupportedDeletions(sv, param=param, pview=pview)
  sv <- removeHemizygous(sv)
  sv <- reviseEachJunction(sv, pview, aview, param)

  sv.revise <- revise(sv, aview, pview, param)
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
  copynumber(sv) <- granges_copynumber(variant(sv), pview)
  calls(sv) <- rpSupportedDeletions(sv, param=param, pview=pview)
  indexImproper(sv) <- updateImproperIndex(sv, maxgap=500)
  calls(sv) <- rpSupportedDeletions(sv, param, pview=pview)
  sv2 <- leftHemizygousHomolog(sv, pview, param)
  sv3 <- rightHemizygousHomolog(sv2, pview, param)
  calls(sv3) <- rpSupportedDeletions(sv3, param, pview)
  message("Refining homozygous boundaries by spanning hemizygous+")
  sv5 <- refineHomozygousBoundaryByHemizygousPlus(sv3)
  sv6 <- callOverlappingHemizygous(sv5)
  sv7 <- removeSameStateOverlapping2(sv6)
  sv <- removeHemizygous(sv7)
  expect_identical(sv, sv.revise)

  i <- order(variant(sv))
  j <- order(variant(expected.sv))
  g <- granges(variant(sv[i]))
  g.expected <- granges(variant(expected.sv[j]))
  expect_identical(start(g), start(g.expected))
  expect_identical(end(g), end(g.expected))
  expect_identical(calls(sv)[i], calls(expected.sv)[j])

  sv.finalize <- finalize_deletions(sv, gr_filters,
                                    pview, bview,
                                    param)
  expect_identical(variant(sv.finalize),
                   variant(sv))

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
