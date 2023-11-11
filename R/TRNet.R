# Purpose: TRNet S4 class of transcriptional regulatory network
# Author: Jielin Yang
# Date: 2023-11-07
# Version: 1.0
# Bugs and Issues: None


#' @title TRNet S4 class of transcriptional regulatory network
#'
#' @keywords internal
#'
#' @description The TRNet S4 class is a class for storing transcriptional
#'              regulatory network. It internally utilizes the igraph package
#'              to store the network, with additional parameters to store
#'              information about the network.
#'
#' @slot network An igraph object containing the transcriptional regulatory
#'               network. The nodes of the network are the genes and the edges
#'               are the regulatory relationships between the genes. This is a
#'               named, directed, and weighted graph.
#' @slot TRNmetadata A data frame containing all interaction data for the
#'                   transcriptional regulatory network in adjacency list
#'                   format. The columns minimally includes "regulator",
#'                   "target", and "regulatorType".
#' @slot predicted A logical value indicating whether the transcriptional
#'                 regulatory network contains inffered interactions.
#' @slot omics A character string indicating the omics data used to construct
#'             the transcriptional regulatory network.
#'
#' @exportClass TRNet
#'
#' @importFrom methods setClass
#'
methods::setClass("TRNet",
  contains = "list",
  slots = c(
    network = "ANY",
    TRNmetadata = "data.frame",
    predicted = "logical",
    omics = "character"
  )
)


#' @title TRNet constructor
#'
#' @description This function is a constructor of the TRNet object
#'
#' @param TRNmetadata A data frame containing all interaction data for the
#'                    transcriptional regulatory network in adjacency list
#'                    format. The columns minimally includes "regulator",
#'                    "target", and "regulatorType".
#' @param predicted A logical value indicating whether the transcriptional
#'                  regulatory network contains inffered interactions.
#' @param omics A character string indicating the omics data used to construct
#'              the transcriptional regulatory network.
#'
#' @return A TRNet object
#'
TRNet <- function(TRNmetadata, predicted, omics) {
  # Check if TRNmetadata is a data frame
  if (!is.data.frame(TRNmetadata)) {
    stop("TRNmetadata must be a data frame")
  } else {
    # Do nothing
  }

  # Check if TRNmetadata contains the required columns
  if (!all(NETWORK_FIELD %in% colnames(TRNmetadata))) {
    stop("TRNmetadata must contain all required fields")
  } else {
    # Do nothing
  }

  # Create TRNet object
  trn <- new("TRNet",
    TRNmetadata = TRNmetadata,
    predicted = predicted,
    omics = omics
  )

  # Create igraph object
  trn <- generatePlot(trn)

  return(trn)
}


#' @title Generate igraph object from TRNet network metadata
#'
#' @rdname TRNet-class
#'
#' @description This function generates an igraph object from the TRNet network
#'              metadata
#'
#' @param trn A TRNet object
#'
#' @return A TRNet object with the igraph object in the network slot
#' 
#' @export
#'
#' @importFrom methods setGeneric setMethod
#' @importFrom igraph graph_from_data_frame
#' @examples
#' # Define some example edges
#' edges <- data.frame(
#'   regulator = c("A", "B", "C"),
#'   target = c("D", "E", "F"),
#'   regulatorType = c("miRNA", "TF", "TF")
#' )
#'
#' # Create TRNet object
#' trn <- TRNet(edges, FALSE, "RNA-seq")
#'
#' # Generate igraph object
#' generatePlot(trn)
#'
methods::setGeneric(
  "generatePlot",
  function(trn) {
    standardGeneric("generatePlot")
  }
)
methods::setMethod(
  "generatePlot", "TRNet",
  function(trn) {
    # Parse vertex metadata
    vertexMetadata <- parseVertexMetadata(trn)
    # Generate igraph object
    network <- igraph::graph_from_data_frame(trn@TRNmetadata[, 1:2],
      directed = TRUE,
      vertices = vertexMetadata
    )

    trn@network <- network
    return(trn)
  }
)


#' @title Parse vertex metadata from TRNet network metadata
#'
#' @rdname TRNet-class
#'
#' @description This function parses vertex metadata from the TRNet network
#'              metadata
#'
#' @param trn A TRNet object
#'
#' @return A data frame containing the vertex metadata
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export
#'
#' @examples
#' # Define some example edges
#' edges <- data.frame(
#'   regulator = c("A", "B", "C"),
#'   target = c("D", "E", "F"),
#'   regulatorType = c("miRNA", "TF", "TF")
#' )
#'
#' # Create TRNet object
#' trn <- TRNet(edges, FALSE, "RNA-seq")
#'
#' # Parse vertex metadata
#' parseVertexMetadata(trn)
#'
methods::setGeneric(
  "parseVertexMetadata",
  function(trn) {
    standardGeneric("parseVertexMetadata")
  }
)
methods::setMethod(
  "parseVertexMetadata", "TRNet",
  function(trn) {
    edgeMetadata <- trn@TRNmetadata
    # Target vertices
    targetVertices <- data.frame(
      name = unique(edgeMetadata$target),
      type = "gene"
    )
    # Regulator vertices
    regulatorVertices <- edgeMetadata[, c("regulator", "regulatorType")]
    colnames(regulatorVertices) <- c("name", "type")
    regulatorVertices <- regulatorVertices[!duplicated(regulatorVertices), ]
    # Combine target and regulator vertices
    vertexMetadata <- rbind(targetVertices, regulatorVertices)
    return(vertexMetadata)
  }
)


#' @title Plot TRNet object
#'
#' @rdname TRNet-class
#'
#' @description This function plots the TRNet object
#'
#' @param trn A TRNet object
#'
#' @return NULL
#'
#' @importFrom methods setGeneric setMethod
#' @importFrom igraph plot.igraph
#'
#' @export
#'
#' @examples
#' # Define some example edges
#' edges <- data.frame(
#'   regulator = c("A", "B", "C"),
#'   target = c("D", "E", "F"),
#'   regulatorType = c("miRNA", "TF", "TF")
#' )
#'
#' # Create TRNet object
#' trn <- TRNet(edges, FALSE, "RNA-seq")
#'
#' # Plot TRNet object
#' plotNetwork(trn)
#'
methods::setGeneric(
  "plotNetwork",
  function(trn) {
    standardGeneric("plotNetwork")
  }
)
methods::setMethod(
  "plotNetwork", "TRNet",
  function(trn) {
    igraph::plot.igraph(trn@network)

  }
)


#' @title Export an igraph object from TRNet object
#' 
#' @rdname TRNet-class
#'
#' @description This function exports an igraph object from TRNet object
#'
#' @param trn A TRNet object
#'
#' @return An igraph object
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export
#'
#' @examples
#' # Define some example edges
#' edges <- data.frame(
#'   regulator = c("A", "B", "C"),
#'   target = c("D", "E", "F"),
#'   regulatorType = c("miRNA", "TF", "TF")
#' )
#'
#' # Create TRNet object
#' trn <- TRNet(edges, FALSE, "RNA-seq")
#'
#' # Export igraph object
#' exportIgraph(trn)
#'
methods::setGeneric(
  "exportIgraph",
  function(trn) {
    standardGeneric("exportIgraph")
  }
)
methods::setMethod(
  "exportIgraph", "TRNet",
  function(trn) {
    return(trn@network)
  }
)


#' @title Export all network interactions from TRNet object
#' 
#' @rdname TRNet-class
#' 
#' @description This function exports all network interactions from TRNet object
#'  
#' @param trn A TRNet object
#' 
#' @return A data frame containing all network interactions
#' 
#' @importFrom methods setGeneric setMethod
#' 
#' @export
#' 
#' @examples
#' # Define some example edges
#' edges <- data.frame(
#'  regulator = c("A", "B", "C"),
#' target = c("D", "E", "F"),
#' regulatorType = c("miRNA", "TF", "TF")
#' )
#' 
#' # Create TRNet object
#' trn <- TRNet(edges, FALSE, "RNA-seq")
#' 
#' # Export all network interactions
#' exportEdgeSet(trn)
#' 
methods::setGeneric(
  "exportEdgeSet",
  function(trn) {
    standardGeneric("exportEdgeSet")
  }
)
methods::setMethod(
  "exportEdgeSet", "TRNet",
  function(trn) {
    return(trn@TRNmetadata)
  }
)


#' @title show method for TRNet object
#' 
#' @rdname TRNet-class
#' 
#' @description This function prints the TRNet object
#' 
#' @param object A TRNet object
#' 
#' @return NULL
#' 
#' @importFrom methods setMethod
#' 
#' @export
#' 
#' @examples
#' # Define some example edges
#' edges <- data.frame(
#'  regulator = c("A", "B", "C"),
#' target = c("D", "E", "F"),
#' regulatorType = c("miRNA", "TF", "TF")
#' )
#' 
#' # Create TRNet object
#' trn <- TRNet(edges, FALSE, "RNA-seq")
#' 
#' # Print TRNet object
#' trn
#' 
methods::setMethod(
  "show", "TRNet",
  function(object) {
    nEdges <- igraph::ecount(object@network)
    nVertices <- igraph::vcount(object@network)
    cat("TRNet object with", nEdges, "edges and", nVertices, "vertices\n")
    nTarget <- length(unique(object@TRNmetadata$target))
    nRegulator <- length(unique(object@TRNmetadata$regulator))
    regTypes <- unique(object@TRNmetadata$regulatorType)
    cat("The network contains", nTarget, "target genes and", nRegulator,
      "regulators\n")
    cat("Including", length(regTypes), "types of regulators:\n")
    for (regType in regTypes) {
      # Find the number of regulatory interactions of each type
      nRegType <- sum(object@TRNmetadata$regulatorType == regType)
      cat("  ", regType, ":", nRegType, "regulatory interactions\n")
    }
    # Predicted interactions
    if (object@predicted) {
      cat("The network is strengthed by co-expression-based predictions\n")
    } else {
      # Do nothing
    }
    # Purely predicted interactions
    predTypes <- setdiff(regTypes, c("miRNA", "TF"))
    if (length(predTypes) > 0) {
      cat("Note: The following types of interactions are purely predicted:\n")
      cat("  ", paste(predTypes, collapse = ", "), "\n")
    } else {
      # Do nothing
    }
    # Omics data
    cat("Omics data involved during network construction:", object@omics, "\n")
    cat("\n")
    cat("A snapshort of network interactions:\n")
    print(head(object@TRNmetadata, 5))
    cat("\n")
    cat("To access all network interactions, use exportEdgeSet()\n")
  }
)


# [END]
