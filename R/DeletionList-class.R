#' @include AllGenerics.R
NULL

setClass("DeletionList", contains="CompressedGRangesList")

#' @rdname DeletionList-methods
#' @aliases DeletionList,list-method
setMethod("DeletionList", "list", function(object){
  DeletionList(GRangesList(object))
})

#' @rdname DeletionList-methods
#' @aliases DeletionList,GRangesList-method
setMethod("DeletionList", "GRangesList", function(object){
  as(object, "DeletionList")
})

#' @rdname DeletionList-methods
#' @aliases DeletionList,missing-method
setMethod("DeletionList", "missing", function(object){
  as(GRangesList(), "DeletionList")
})


#' @rdname splitReads
#' @aliases splitReads,DeletionList,GRangesList-method
setReplaceMethod("splitReads", c("DeletionList", "GRangesList"),
                 function(object, value){
                   orig_order <- names(object)
                   object2 <- object[ names(object) %in% names(value) ]
                   object2 <- object2 [ names(value) ]
                   for (i in seq_len(length(object2))) {
                     splitReads(object2[[i]]) <- value[[i]]
                   }
                   notchanged <- object [ !names(object) %in% names(object2) ]
                   if(length(notchanged) > 0){
                     object3 <- c(notchanged, object2)
                   } else object3 <- object2
                   object3 <- object3[ orig_order ]
                   object3
                 })


#' @rdname splitReads
#' @aliases splitReads,DeletionList-method
setMethod("splitReads", "DeletionList", 
          function(object){
            split_reads <- vector("list", length(object))
            for(i in seq_along(object)){
              split_reads[[i]] <- splitReads(object[[i]])
            }
            grl <- GRangesList(split_reads)
            names(grl) <- names(object)
            grl
          })
