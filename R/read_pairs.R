## a convenience function for unit tests. Not useful in practice because these
## functions may take a long time to run on a large BAM file
##
## 
listReadPairs <- function(bam.file, segments, build){
  dparam <- DeletionParam()
  aparam <- ampliconParams()
  is.amp <- segments$seg.mean > aparam$AMP_THR
  is.del <- segments$seg.mean < hemizygousThr(dparam)
  amp.gr <- segments[is.amp]
  del.gr <- segments[is.del]
  proper.amp <- get_readpairs2(reduce(amp.gr,
                                      min.gapwidth=2000),
                               bam.file)
  proper.del <- properReadPairs(bam.file,
                                gr=reduce(del.gr, min.gapwidth=2000),
                                param=dparam)


  irp.params <- improperAlignmentParams(mapqFilter=30)
  improper_rp <- getImproperAlignmentPairs(bam.file, irp.params, build=build)
  read_pairs <- list(proper_del=proper.del,
                     proper_amp=proper.amp,
                     improper=improper_rp)
}
