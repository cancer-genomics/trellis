context("vignette")

## Create a unit test that will fail
test_that("vignette deletion", {
    tmp  <- readRDS("temp.rds")
    names(tmp)
    pdata <- tmp$pdata
    dp <- tmp$dp
    deletions <- sv_deletions(preprocess=pdata, param=dp)
    L <- length(variant(deletions))
    expect_identical(L, 4L)
})


test_that("step sv_deletion", {
    tmp  <- readRDS("temp.rds")
    preprocess <- tmp$pdata
    param <- tmp$dp
    gr_filters <- genomeFilters(preprocess$genome)
    segs <- preprocess$segments
    preprocess$segments <- segs[segs$seg.mean < hemizygousThr(param)]
    sv <- deletion_call(preprocess, gr_filters, param)
    expect_identical(length(sv), 4L)
    expect_equal(variant(sv)$seg.mean[4], -8.6, tolerance=0.05)
    calls(sv) <- rpSupportedDeletions(sv, param=param, bins=preprocess$bins)
    expect_identical(calls(sv)[4], "homozygous+")    
    ##
    ## Define region of interest (ROI)
    ##
    ROI <- GRanges("chr15", IRanges(63201001, 63209262),
                   seg.mean=-8.6549, sample="CGOV44T")
    seqlevels(ROI) <- seqlevels(variant(sv))
    seqinfo(ROI) <- seqinfo(variant(sv))
    names(ROI) <- "sv4"
    
    sv <- removeHemizygous(sv)
    ## Extracting improper read pairs with a high MAPQ
    improper_rp <- preprocess$read_pairs[["improper"]]
    ## avoid namespace issues with dplyr
    first <- GenomicAlignments::first
    last <- GenomicAlignments::last
    F <- first(improper_rp); L <- last(improper_rp)
    mapq <- mcols(F)$mapq > 30 & mcols(L)$mapq > 30
    improper_rp <- improper_rp[mapq]
    ##
    ## Check that there are improper read pairs within 1kb of the ROI
    ##
    irp_near_deletion <- subsetByOverlaps(improper_rp, ROI, maxgap=1e3)
    expect_identical(length(irp_near_deletion), 55L)
    ## Revise deletion boundaries using improper read pairs
    ## and proper read pairs
    if(FALSE){
        rej_args <- list(sv=sv, improper_rp=improper_rp)
        saveRDS(rej_args, file="temp_reviseEachJunction_args.rds")
    }
    sv <- reviseEachJunction(sv=sv,
                             bins=preprocess$bins,
                             improper_rp=improper_rp,
                             param=param)
    ##Why are improper read-pairs not added to the SV object here?
    expect_true(length(improper(sv)) > 0)    
    ##
    ## Next, step through .reviseEachJunction
    ##
})

test_that("step .reviseEachJunction", {
    tmp  <- readRDS("temp.rds")
    preprocess <- tmp$pdata
    bins <- preprocess$bins
    param <- tmp$dp    
    rej_args  <- readRDS("temp_reviseEachJunction_args.rds")
    sv <- rej_args$sv
    improper_rp <- rej_args$improper_rp
    svs <- .reviseEachJunction(sv, bins, improper_rp, param)
    expect_true(length(improper(svs)) > 0)
    ## fails even when improper read pairs are manually added
    sv2 <- sv
    improper(sv2) <- improper_rp
    ## note, the improper<- assignment method only adds improper read pairs near the candidate deletions
    ## -- this is reasonable
    length(improper(sv2)) < length(improper_rp)
    irp_near_deletion <- subsetByOverlaps(improper(sv2), ROI, maxgap=1e3)
    expect_identical(length(irp_near_deletion), 55L)
    ##
    ## Below, we dig into .reviseEachJunction to understand
    ## why the improper read pairs are dropped
    ##
    svs <- reviseDeletionBorders(sv)
    ##
    ## for the duplicated ranges, revert back to the original
    ##
    svs[duplicated(svs)] <- variant(sv)[duplicated(svs)]
    variant(sv) <- svs
    copynumber(sv) <- granges_copynumber2(svs, bins)
    ## Revise homozygous deletion borders using
    ## the absence of proper read pairs (this likely only works)
    ## with 100% pure samples (e.g. cell lines)
    variant(sv) <- homozygousBorders(sv, svs)
    is.dup <- duplicated(variant(sv))
    if(any(is.dup)){
      sv <- sv[!is.dup]
    }
    ## For hemizygous deletions, look for gaps in proper read pairs. If there
    ## is a gap, create a homozygous deletion in the gap and add to the StructuralVariant object.
    sv <- hemizygousBorders(sv, param)
    ## Update the improper read pairs supporting deletions 
    ## since deletion boundaries have been updated
    irp <- improperRP(variant(sv), improper_rp, param=param)
    improper(sv) <- irp
    expect_equal(length(improper(sv)), 74L)
    if(FALSE){
        tmp <- list(gr=variant(sv), irp=irp)
        saveRDS(tmp, file="temp_updateImproperIndex2_args.rds")
    }
    indexImproper(sv) <- updateImproperIndex2(variant(sv), irp, maxgap=500)
    ##
    ## Indexing the improper read pairs ends up dropping the reads
    ## spanning the homozygous deletion.
    ##
    expect_equal(length(improper(sv)), 51L) 
    ##
    ## Next, step through updateImproperIndex2 function
    ##
})

test_that("step updateImproperIndex2", {
    tmp <- readRDS("temp_updateImproperIndex2_args.rds")
    gr <- tmp$gr
    irp <- tmp$irp
    maxgap <- 2e3
    left_boundary <- resize(gr, width=2)
    right_boundary <- resize(gr, width=2, fix="end")
    index_all <- setNames(vector("list", length(gr)), names(gr))
    if(FALSE){
        ## simpler??
        for(i in seq_along(gr)){
            tmp <- (overlapsAny(first(irp), left_boundary[i], maxgap=maxgap) &
                    overlapsAny(last(irp), right_boundary[i], maxgap=maxgap)) |
                (overlapsAny(first(irp), right_boundary[i], maxgap=maxgap) &
                 overlapsAny(last(irp), left_boundary[i], maxgap=maxgap))
            index_all[[i]] <- which(tmp)
        }
        index_irp <- index_all
        return(index_irp)
    }
    hitsLeft <- findOverlaps(left_boundary, irp, maxgap=maxgap)
    hitsRight <- findOverlaps(right_boundary, irp, maxgap=maxgap)

    index_left <- split(subjectHits(hitsLeft),
                        names(left_boundary)[queryHits(hitsLeft)])
    index_right <- split(subjectHits(hitsRight),
                         names(right_boundary)[queryHits(hitsRight)])

    ## concatenate indices for each element of the index list
    index_all <- setNames(vector("list", length(gr)), names(gr))
    index_right2 <- index_left2 <- index_all
    i <- match(names(index_left), names(index_all))
    index_left2[i] <- index_left
    i <- match(names(index_right), names(index_all))
    index_right2[i] <- index_right
    updated_index_list <- mapply(function(x,y) unique(intersect(x,y)),
                                 index_left2, index_right2)
    ##
    ## some of the elements in this list may be NULL
    ##
    ## Assess whether there are any rearranged read pair clusters
    ## 1. order read pairs by left most alignment
    irp_id <- lapply(updated_index_list, function(i, irp) names(irp)[i], irp=irp)
    lrp <- leftAlwaysFirst(irp) ##, names(irp)
    names(lrp) <- names(irp)
    ##
    ## lrp is no longer named because of recent change to firstIsLeft
    ##
    expect_identical(names(irp), names(lrp))
    lrplist <- lapply(irp_id, function(id, lrp) lrp[names(lrp) %in% id], lrp=lrp)
    expect_true(length(lrplist[[2]]) > 0)
    ## 2. cluster the read pairs for each element
    cl_list <- lapply(lrplist, clusterReadPairs)
    ## all read pairs belong to one cluster
    expect_true(all(cl_list[[1]] == 0))
    ## 3. identify id of cluster with most read pairs
    clid <- lapply(cl_list, function(cl) names(which.max(table(cl))))
    ## 4. extract the names of the improper read pairs that correspond
    ## to the above cluster
    irp_id2 <- mapply(function(lrp, clid, cl) names(lrp)[cl==clid],
                      lrp=lrplist, clid=clid, cl=cl_list)
    ## Update the index a second time to include only those improper
    ## read pairs belonging to the cluster
    index_irp <- lapply(irp_id2, function(id, irp) match(id, names(irp)), irp=irp)
    ## the NULLs are now converted to NAs.
    ## Subsetting a vector by NULL returns an empty vector. Convert NAs back to nulls
    na_index <- which(sapply(index_irp, function(x) any(is.na(x))))
    for(j in na_index){
      index_irp[j] <- list(NULL)
    }
    expect_identical(length(index_irp[[2]]), 51L)
})
