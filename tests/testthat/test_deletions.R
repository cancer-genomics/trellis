context("Call deletions")
test_that("deletions", {
  library(svovarian)
  library(svfilters)
  library(svcnvs)
  id <- "CGOV21T"
  dp <- projectOvarian(rootname="OvarianData2")
  data(lymphoblast_filters_hg19)
  data(lowMappabilityBins_hg19)
  data(binAssemblyGaps_hg19)
  gfilters <- as(lymphoblast_filters_hg19, "list")
  gfilters$map <- lowMappabilityBins_hg19
  gfilters$gc <- binAssemblyGaps_hg19
  germline_filters <- reduce(unlist(GRangesList(lapply(gfilters, granges))))
  bviews <- readRDS(file.path(dp[1], "bviews_hg19.rds"))
  del_file <- file.path(dp["1deletions"], rdsId(bviews[, id]))
  file2 <- gsub(".rds", "_copy.rds", del_file)
  if(file.exists(del_file)){
    unlink(file2)
    file.copy(del_file, file2)
    unlink(del_file)
  }
  grl <- readRDS(file.path(dp["segment"], "grl_hg19.rds"))
  gr <- grl[[id]]
  gr <- gr[gr$seg.mean < log2(0.75)]
  grl_del <- GRangesList(list(CGOV21T=gr))
  sv_dels <- sv_deletion_exp(dirs=dp, grl=grl_del[id],
                             bviews=bviews[, id],
                             gr_filters=germline_filters)
  expect_is(sv_dels[[1]], "StructuralVariant")
  expect_true(file.exists(del_file))
  unlink(del_file)
  file.copy(file2, del_file)
  unlink(file2)
})


test_that("deletions2", {
  library(svovarian)
  library(Rsamtools)
  library(svpreprocess)
  library(svcnvs)
  dirs <- projectOvarian()
  ut_path <- file.path(dirs[["unit_test"]], "deletion_pipeline")
   resaveObject2 <- function(path){
     object <- readRDS(path)
     cl <- attributes(object)$class
     attributes(cl)$package <- "svclasses"
     attributes(object)$class <- cl
     saveRDS(object, file=path)
   }
  del_list <- readRDS(file.path(ut_path, "del_list.rds"))
  bviews <- readRDS(file.path(dirs[1], "bviews_hg19.rds"))
  index <- match(gsub(".bam", "", names(del_list)[12]), bamSamples(bviews)$Flow.Cell.ID)
  bview <- bviews[, index]
  ## "CGOV12T"
  cnv_start <- del_list[[12]]
  ## after some filtering, we have an updated object
  saved_result <- readRDS(file.path(ut_path, "cnv1.rds"))

  ##
  ## Reproduce saved results
  ##
  data(all_filters_hg19, package="Ovarian")
  param <- DeletionParam()
  germline_filters <- all_filters
  pview <- sv_preprocess(bview, dirs)
  names(cnv_start) <- paste0("sv", seq_along(cnv_start))
  test <- germlineFilters(cnv_start, all_filters, pview)
  expect_true(all(overlapsAny(saved_result, cnv_start)))
  expect_true(all(overlapsAny(saved_result, test)))
  expect_true(all(overlapsAny(test, saved_result)))
  ##
  ##
  ##
  aview <- AlignmentViews2(bview, dirs)  
  REMOTE <- file.exists(bamPaths(bview))
  if(REMOTE){
    ## On cluster
    sv_saved <- readRDS(file.path(ut_path, "sv1.rds"))
    sv_test <- .deletion_call(aview, pview, cnv_start, all_filters)
    expect_identical(length(proper(sv_test)), length(proper(sv_saved)))
    expect_identical(length(improper(sv_test)), length(improper(sv_saved)))
    expect_identical(calls(sv_test), calls(sv_saved))
  }
  sv_saved <- readRDS(file.path(ut_path, "sv2.rds"))
  ## Reproduce initial calls
  sv_input <- readRDS(file.path(ut_path, "sv1.rds"))
  cncalls <- rpSupportedDeletions(sv_input, param=param, pview=pview)  
  calls(sv_input) <- cncalls
  expect_identical(sv_saved, sv_input)

  ## Result after revising junctions
  if(REMOTE){
    is_hemizygous <- calls(sv_input)=="hemizygous"
    sv <- sv_input[!is_hemizygous]
    test <- reviseEachJunction(sv, pview, aview, param)
    answer <- readRDS(file.path(ut_path, "sv3.rds"))
    expect_identical(length(proper(test)), length(proper(answer)))
    expect_identical(length(improper(test)), length(improper(answer)))
    expect_identical(calls(test), calls(answer))
  }

  sv_test <- readRDS(file.path(ut_path, "sv3.rds"))
  copynumber(sv_test) <- granges_copynumber(variant(sv_test), pview)
  saved_result <- readRDS(file.path(ut_path, "sv4.rds"))
  expect_identical(copynumber(sv_test), copynumber(saved_result))
  
  calls(sv_test) <- rpSupportedDeletions(sv_test, param=param, pview=pview)  
  saved_result <- readRDS(file.path(ut_path, "sv4.rds"))
  expect_identical(calls(sv_test), calls(saved_result))

  is_hemizygous <- calls (sv_test) == "hemizygous"
  sv_test <- sv_test[!is_hemizygous]
  sv_test <- removeSameStateOverlapping(sv_test)
  saved_result <- readRDS(file.path(ut_path, "sv5.rds"))
  expect_equal(sv_test, saved_result)


  indexImproper (sv_test) <- updateImproperIndex (sv_test, maxgap=500)
  calls (sv_test) <- rpSupportedDeletions (sv_test, param)
  is_hemizygous <- calls (sv_test) == "hemizygous"
  sv_test <- sv_test[!is_hemizygous]
  sv_test2 <- SVFilters(sv_test, all_filters, pview, param=param)
  saved_result <- readRDS(file.path(ut_path, "sv6.rds"))
  expect_equal(sv_test2, saved_result)

  if(REMOTE){
    sv_test <- groupSVs(sv_test2)
    id <- names(aview)
    sv_test <- allProperReadPairs(sv_test, param, bfile=bamPaths(bview), zoom.out=1)
    sv_result <- readRDS(file=file.path(ut_path, "sv7.rds"))
    expect_equal(sv_test, sv_result)
  }
  if(REMOTE){
    sv_test <- sv_deletions(cnv_start, aview, bview,
                            pview,
                            gr_filters=all_filters,
                            param=DeletionParam())
    expect_equivalent(sv_test, sv_result)
  }
  ## because the above test on REMOTE works, we only need to verify
  ## that sv7.rds is in fact the list of tabled deletions read in by
  ## the function readDeletions.
  sv_result <- readRDS(file=file.path(ut_path, "sv7.rds"))
  dels <- readDeletions(dirs, seqinfo(pview))
  test_cgov12 <- variant(sv_result)
  names(test_cgov12) <- NULL
  saved_dels <- dels[dels$id == "CGOV12T"]
  expect_identical(length(saved_dels), length(test_cgov12))
  expect_identical(granges(saved_dels),
                   granges(test_cgov12))
})

test_that("sv_deletion_exp", {
  library(svovarian)
  library(Rsamtools)
  library(svpreprocess)
  dirs <- projectOvarian()
  ut_path <- file.path(dirs[["unit_test"]], "deletion_pipeline")
  del_list <- readRDS(file.path(ut_path, "del_list.rds"))
  bviews <- readRDS(file.path(dirs[1], "bviews_hg19.rds"))
  pviews <- sv_preprocess(bviews, dirs)
  index <- match(gsub(".bam", "", names (del_list)[12]),
                 bamSamples (bviews)$Flow.Cell.ID)
  bview <- bviews[, index]
  aview <- AlignmentViews2(bview, dirs)  
  REMOTE <- file.exists(bamPaths(bview))
  if(REMOTE){
    ## "CGOV12T"
    del_list <- del_list[12]
    names(del_list) <- "CGOV12T"
    data(all_filters_hg19, package="Ovarian", envir=environment())
    dirs[["deletions"]] <- tempdir()
    test <- sv_deletion_exp(dirs=dirs,
                            grl=del_list,
                            bviews=bview,
                            gr_filters=all_filters)
    test <- test[["CGOV12T"]]
    sv7 <- readRDS(file=file.path(ut_path, "sv7.rds"))
    expect_equivalent(test, sv7)
  }
})

test_that("unlist_sv_dels", {
  library(svovarian)
  dirs <- projectOvarian()
  bviews <- readRDS(file.path(dirs[1], "bviews_hg19.rds"))
  saved_gr <- readDeletions(dirs, seq_info=seqinfo(bamRanges(bviews)))
  del_list <- readRDS(file.path(dirs[["unit_test"]], "deletion_list.rds"))
  bviews <- bviews[, unique(saved_gr$id)]
  sv_dels <- sv_deletion_exp(dirs=dirs, bviews=bviews)
  gr_dels <- unlist_sv_dels(sv_dels)
  expect_is(gr_dels, "GRanges")
})

test_that("Compare to saved deletions", {
  library(svovarian)
  dirs <- projectOvarian()
  bviews <- readRDS(file.path(dirs[1], "bviews_hg19.rds"))
  saved_gr <- readDeletions(dirs, seq_info=seqinfo(bamRanges(bviews)))
  del_list <- readRDS(file.path(dirs[["unit_test"]], "deletion_list.rds"))
  
  bviews <- bviews[, unique(saved_gr$id)]
  sv_dels <- sv_deletion_exp(dirs=dirs, bviews=bviews)
  del_gr <- unlist_sv_dels(sv_dels)

  ## appear to be a bit different
  ##  in_common <- intersect(del_gr, saved_gr)
  ##  frac <- sum(width(in_common))/sum(width(del_gr)) ## 0.76
  ##  frac <- sum(width(in_common))/sum(width(saved_gr)) ## 0.76
  sum(width(del_gr))/1e6 ## 207Mb
  sum(width(saved_gr))/1e6 ## 204Mb

  mn <- mean(overlapsAny(del_gr, saved_gr))
  expect_more_than(mn, 0.98)
  mn <- mean(overlapsAny(saved_gr, del_gr))
  expect_more_than(mn, 0.995)

  saved_amps <- readAmps(seqinfo=seqinfo(bamRanges(bviews)))
  ag_list <- sv_amplicon_exp(dirs, bviews)
  amp_gr <- unlist_sv_amps(ag_list)
  
  mn <- mean(overlapsAny(amp_gr, saved_amps))  ## this does not verify concordance of the same samples...
  expect_more_than(mn, 0.99)
  mn <- mean(overlapsAny(saved_amps, amp_gr))
  expect_more_than(mn, 0.995)
})

