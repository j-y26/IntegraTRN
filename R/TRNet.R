# Purpose: TRNet S4 class of transcriptional regulatory network
# Author: Jielin Yang
# Date: 2023-11-07
# Version: 1.0
# Bugs and Issues: None


#' @title TRNet S4 class of transcriptional regulatory network
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
#' @slot .Data A list containing additional network information
#'
#' @exportClass TRNet
#'
setClass("TRNet",
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
#' @references
#' \insertRef{csardi2006igraph}{IntegraTRN}
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
#' @importFrom igraph graph_from_data_frame
#'
#' @references
#' \insertRef{csardi2006igraph}{IntegraTRN}
#'
setGeneric(
  "generatePlot",
  function(trn) {
    standardGeneric("generatePlot")
  }
)
setMethod(
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
#' @aliases parseVertexMetadata,TRNet-method
#'
#' @description This function parses vertex metadata from the TRNet network
#'              metadata
#'
#' @param trn A TRNet object
#'
#' @return A data frame containing the vertex metadata
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
setGeneric(
  "parseVertexMetadata",
  function(trn) {
    standardGeneric("parseVertexMetadata")
  }
)
setMethod(
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
#' @aliases plotNetwork,TRNet-method
#'
#' @description This function plots the TRNet object
#'
#' @param trn A TRNet object
#' @param interactive A logical value indicating whether the plot should be
#'                    interactive. If TRUE, an interactive plot will be
#'                    plotted using the networkD3 package. If FALSE, a static
#'                    plot will be plotted using the igraph package.
#' @param ... Additional arguments passed to the plot function for plotting
#'            a static plot using the igraph package, or the forceNetwork
#'            function for interactive plots using the networkD3. See the
#'            respective package documentation for details.
#'
#' @return NULL
#' @importFrom igraph plot.igraph layout_with_drl categorical_pal
#' @importFrom networkD3 igraph_to_networkD3 forceNetwork
#'
#' @export
#'
#' @references
#' \insertRef{csardi2006igraph}{IntegraTRN}
#'
#' \insertRef{networkd3}{IntegraTRN}
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
#' # Plot interactive TRNet object
#' plotNetwork(trn, interactive = TRUE)
#'
setGeneric(
  "plotNetwork",
  function(trn, interactive = FALSE, ...) {
    standardGeneric("plotNetwork")
  }
)
setMethod(
  "plotNetwork", "TRNet",
  function(trn, interactive = FALSE, ...) {
    network <- trn@network
    # Label vertices with their type
    vertexMetadata <- parseVertexMetadata(trn)
    vertexMetadata <- stats::setNames(vertexMetadata$type, vertexMetadata$name)
    # Plot network
    if (interactive == TRUE) {
      # Convert igraph object to networkd3 object
      network <- networkD3::igraph_to_networkD3(network, group = vertexMetadata)
      networkD3::forceNetwork(
        Links = network$links,
        Nodes = network$nodes,
        Source = "source",
        Target = "target",
        NodeID = "name",
        Group = "group",
        opacity = 0.8,
        zoom = TRUE,
        ...
      )
    } else {
      # Assign vertex color based on type of gene/transcript/TF
      colors <- igraph::categorical_pal(length(unique(vertexMetadata)))
      colors <- stats::setNames(colors, unique(vertexMetadata))
      igraph::V(network)$color <- colors[igraph::V(network)$type]
      # Define layout
      layout <- igraph::layout_with_drl(network)
      # Plot static network
      igraph::plot.igraph(network,
        layout = layout,
        vertex.color = igraph::V(network)$color,
        vertex.label = igraph::V(network)$name,
        vertex.label.color = "black",
        vertex.label.cex = 0.8,
        vertex.label.dist = 0.5,
        edge.arrow.size = 0.3,
        edge.arrow.width = 0.3,
        edge.width = 0.5,
        edge.color = "grey",
        ...
      )
    }
  }
)


#' @title Export an igraph object from TRNet object
#'
#' @aliases exportIgraph,TRNet-method
#'
#' @description This function exports an igraph object from TRNet object
#'
#' @param trn A TRNet object
#'
#' @return An igraph object
#'
#' @export
#'
#' @references
#' \insertRef{csardi2006igraph}{IntegraTRN}
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
setGeneric(
  "exportIgraph",
  function(trn) {
    standardGeneric("exportIgraph")
  }
)
setMethod(
  "exportIgraph", "TRNet",
  function(trn) {
    return(trn@network)
  }
)


#' @title Export all network interactions from TRNet object
#'
#' @aliases exportEdgeSet,TRNet-method
#'
#' @description This function exports all network interactions from TRNet object
#'
#' @param trn A TRNet object
#'
#' @return A data frame containing all network interactions
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
#' # Export all network interactions
#' exportEdgeSet(trn)
#'
setGeneric(
  "exportEdgeSet",
  function(trn) {
    standardGeneric("exportEdgeSet")
  }
)
setMethod(
  "exportEdgeSet", "TRNet",
  function(trn) {
    return(trn@TRNmetadata)
  }
)


#' @title show method for TRNet object
#'
#' @aliases show,TRNet-method
#'
#' @description This function prints the TRNet object
#'
#' @param object A TRNet object
#'
#' @return NULL
#' @importFrom igraph ecount vcount
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
#' # Print TRNet object
#' trn
#'
setMethod(
  "show", "TRNet",
  function(object) {
    nEdges <- igraph::ecount(object@network)
    nVertices <- igraph::vcount(object@network)
    cat("TRNet object with", nEdges, "edges and", nVertices, "vertices\n")
    nTarget <- length(unique(object@TRNmetadata$target))
    nRegulator <- length(unique(object@TRNmetadata$regulator))
    regTypes <- unique(object@TRNmetadata$regulatorType)
    cat(
      "The network contains", nTarget, "target genes and", nRegulator,
      "regulators\n"
    )
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
    print(utils::head(object@TRNmetadata, 5))
    cat("\n")
    cat("To access all network interactions, use exportEdgeSet()\n")
  }
)


# [END]
