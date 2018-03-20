#' @include AllGenerics.R
NULL

#' @rdname AmpliconGraph-class
#' @export
#' @param amplicons genomic intervals of amplicons as a \code{GRanges} object
ampliconNames <- function(amplicons){
  if(length(amplicons)==0) return(character())
  gsub(" ", "", paste(seqnames(amplicons), ":",
                      prettyNum(start(amplicons), big.mark=","), sep=""))
}


#' A class for representing somatic amplicons
#'
#' Somatic amplicons are often grouped by a bridge-breakage-fusion
#' cycles.  We represent the set of amplicons identified for a sample
#' as a graph with amplicons as the nodes. An edge drawn between two
#' nodes indicates that five or more read pairs with aberrant spacing
#' support the link between these amplicons.  As the class is an
#' extension of the \code{graphNEL} class, plot methods to visualize
#' the graph are readily available.
#'
#' @seealso See \code{sv_amplicon_exp} for generating a list of
#'   \code{AmpliconGraphs} for an experiment.
#' @rdname AmpliconGraph-class
#' @export
setClass("AmpliconGraph", representation(graph="graphNEL",
                                         ranges="GRanges",
                                         border_size="numeric",
                                         query="GRanges",
                                         assembly_gaps="GRanges",
                                         germline_cnv="GRanges",
                                         outliers="GRanges"))


setValidity("AmpliconGraph", function(object){
  msg <- TRUE
  if(length(ampliconRanges(object)) > 0 || length(nodes(object)) > 0){
    if(length(ampliconRanges(object)) != length(nodes(object))){
      msg <- "ampliconRanges should be the same length as nodes"
      return(msg)
    }
    if(!identical(names(ampliconRanges(object)), nodes(object))){
      msg <- "names of ampliconRanges should be identical to the nodes labels"
      return(msg)
    }
  }
  if(any(names(ranges(object))=="") || is.null(names(ranges(object)))){
    msg <- "GRanges in slot ranges must be named.  See ampliconNames()"
    return(msg)
  }
  msg
})

totalWidth <- function(object) sum(as.numeric(width(object)))

setMethod("show", "AmpliconGraph", function(object){
  cat("An object of class `AmpliconGraph'\n")
  cat("   number of amplicons (nodes):",
      numNodes(object@graph), "\n")
  cat("   number of edges:", numEdges(object@graph), "\n")
  cat("   number of non-amplicon ranges:", sum(!isAmplicon(object)), "\n")
  sizeofamp <- totalWidth(ampliconRanges(object))
  cat("   size of amplicons:", round(sizeofamp/1e6, 2), "Mb \n")
  cat("   size of non-amplicon ranges:",
      round(totalWidth(nonampliconRanges(object))/1e6, 2), "Mb \n")
  filtersize <- round(totalWidth(assemblyGaps(object))/1e6,2)
  cat("   size of filters (centromes/heterochromatin + assembly gaps):",
      filtersize, "Mb \n")
  sizeofquery <- totalWidth(queryRanges(object))
  cat("   size of query:", round(sizeofquery/1e6, 2), "Mb \n")
})


node2 <- function(name, sep="-") sapply(name, function(x) strsplit(x, sep)[[1]][2])
node1 <- function(name, sep="-") sapply(name, function(x) strsplit(x, sep)[[1]][1])

#' @rdname AmpliconGraph-class
#' @aliases addNode,character,AmpliconGraph-method
#' @param node name of amplicon 
#' @param edges name of edge (e.g., concatenate two amplicon names)
setMethod("addNode", c("character", "AmpliconGraph"),
          function(node, object, edges){
            index <- match(node, names(ranges(object)))
            new_ranges <- ranges(object)[index]
            ranges(object)$is_amplicon[index] <- TRUE
            all_amplicons <- sort(c(ampliconRanges(object),
                                    new_ranges))
            all_amplicons <- all_amplicons[!duplicated(names(all_amplicons))]
            existing_edges <- gsub("~", "-", edgeNames(graph(object)))
            edges <- c(existing_edges, edges)
            graph(object) <- graphNEL(nodes=names(all_amplicons))
            graph(object) <- addEdge(node1(edges),
                                     node2(edges),
                                     graph(object))
            object
          })

#' @rdname AmpliconGraph-class
#' @aliases amplicons,AmpliconGraph-method
setMethod(amplicons, "AmpliconGraph", function(object) {
  ranges(object)[names(ampliconRanges(object)), ]
})

#' @rdname AmpliconGraph-class
#' @aliases amplicons<-,AmpliconGraph,ANY-method
setReplaceMethod("amplicons", "AmpliconGraph", function(object, value){
  object@amplicons <- value
  object
})

#' @rdname AmpliconGraph-class
#' @aliases assemblyGaps,AmpliconGraph-method
setMethod(assemblyGaps, "AmpliconGraph", function(object) object@assembly_gaps)

#' @rdname AmpliconGraph-class
#' @aliases assemblyGaps<-,AmpliconGraph,ANY-method
#' @param value genomic amplicon intervals as a \code{GRanges} object
setReplaceMethod("ampliconRanges", "AmpliconGraph", function(object, value){
  amplicons(object)@ranges <- value
  object
})

#' @rdname AmpliconGraph-class
#' @aliases ampliconRanges,AmpliconGraph-method
setMethod("ampliconRanges", "AmpliconGraph", function(object){
  ar <- ampliconRanges(ranges(object))
  ar[names(ar) %in% nodes(object)]
})

#' @rdname AmpliconGraph-class
#' @aliases ampliconRanges,GRanges-method
setMethod(ampliconRanges, "GRanges", function(object){
  object[isAmplicon(object) & isSomatic(object)]
})

#' @rdname AmpliconGraph-class
#' @aliases graph,AmpliconGraph-class
setMethod(graph, "AmpliconGraph", function(object) object@graph)

#' @rdname AmpliconGraph-class
#' @aliases driver,AmpliconGraph-method
setMethod(driver, "AmpliconGraph", function(object){
  ampliconRanges(object)$driver
})

#' @rdname AmpliconGraph-class
#' @aliases edges,AmpliconGraph-method
 #' @param which ignored
setMethod("edges", "AmpliconGraph", function(object, which, ...){
  edges(graph(object))
})

#' @rdname AmpliconGraph-class
#' @aliases graph,AmpliconGraph,ANY-method
setReplaceMethod("graph", "AmpliconGraph", function(object, value){
  object@graph <- value
  object
})

#' @rdname AmpliconGraph-class
#' @aliases isAmplicon,AmpliconGraph-method
setMethod("isAmplicon", "AmpliconGraph", function(object) {
  isAmplicon(ranges(object)) & names(ranges(object)) %in% nodes(object)
})

#' @rdname AmpliconGraph-class
#' @aliases isAmplicon,GRanges-method
setMethod("isAmplicon", "GRanges", function(object){
  if(length(object)==0) return(logical())
  object$is_amplicon
})

#' @rdname AmpliconGraph-class
#' @aliases isSomatic,GRanges-method
setMethod("isSomatic", "GRanges", function(object){
  if(length(object)==0) return(logical())
  object$overlaps_germline != "fully_germline"
})


#' @export
#' @rdname AmpliconGraph-class
#' @aliases length,AmpliconGraph-method
setMethod("length", "AmpliconGraph", function(x) length(nodes(x)))

#' @rdname AmpliconGraph-class
#' @aliases nodes,AmpliconGraph-method
#' @param ... ignored
setMethod("nodes", "AmpliconGraph", function(object) nodes(graph(object)))
##setMethod("nodes", "AmpliconGraph", function(object, ...)   names(ampliconRanges(object)))

#' @rdname AmpliconGraph-class
#' @aliases names<-,AmpliconGraph,ANY-method
#' @param x a \code{AmpliconGraph}
setReplaceMethod("names", "AmpliconGraph", function(x, value){
  x@id <- value
  x
})

#' @rdname AmpliconGraph-class
#' @aliases nonAmpliconRanges,AmpliconGraph-method
setMethod(nonampliconRanges, "AmpliconGraph", function(object) {
  nonampliconRanges(ranges(object))
})

#' @rdname AmpliconGraph-class
#' @aliases nonAmpliconGranges,GRanges-method
setMethod(nonampliconRanges, "GRanges", function(object){
  object[!object$is_amplicon]
})

#' @rdname AmpliconGraph-class
#' @aliases numEdges,AmpliconGraph-method
setMethod("numEdges", "AmpliconGraph", function(object) numEdges(graph(object)))

#' @rdname AmpliconGraph-class
#' @aliases numNodes,AmpliconGraph-method
setMethod("numNodes", "AmpliconGraph", function(object) numNodes(graph(object)))


#' @rdname AmpliconGraph-class
#' @aliases queryRanges,AmpliconGraph-method
setMethod("queryRanges", "AmpliconGraph", function(object) object@query)

#' @rdname AmpliconGraph-class
#' @aliases queryRanges,AmpliconGraph,ANY-method
setReplaceMethod("queryRanges", "AmpliconGraph", function(object, value) {
  object@query <- value
  object
})

#' @rdname AmpliconGraph-class
#' @aliases ranges,AmpliconGraph-method
setMethod("ranges", "AmpliconGraph", function(x, ...) x@ranges)

#' @rdname AmpliconGraph-class
#' @aliases ranges,AmpliconGraph,ANY-method
setReplaceMethod("ranges", "AmpliconGraph", function(x, value){
  x@ranges <- value
  x
})
