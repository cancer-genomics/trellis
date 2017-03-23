## a convenience function for unit tests
listReadPairs <- function(bam.file, segments){
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
  improper_rp <- getImproperAlignmentPairs(bamfile, irp.params)
  read_pairs <- list(proper_del=proper.del,
                     proper_amp=proper.amp,
                     improper=improper_rp)
}
