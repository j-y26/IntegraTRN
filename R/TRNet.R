# Purpose: TRNet S4 class of transcriptional regulatory network
# Author: Jielin Yang
# Date: 2023-11-07
# Version: 1.0
# Bugs and Issues: None


# Register S3 igraph class
methods::setOldClass("igraph")


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
#' @importFrom methods setClass setOldClass
#' @importClassesFrom igraph graph
#' 
methods::setClass("TRNet",
                  contains = "list",
                  slots = c(network = "igraph",
                            TRNmetadata = "data.frame",
                            predicted = "logical",
                            omics = "character"))


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
  if (!all(NETOWRK_FIELD %in% colnames(TRNmetadata))) {
    stop("TRNmetadata must contain all required fields")
  } else {
    # Do nothing
  }
  
  # Create TRNet object
  trn <- new("TRNet",
               TRNmetadata = TRNmetadata,
               predicted = predicted,
               omics = omics)

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
#' @importFrom igraph graph_from_data_frame
#' 
#' @examples
#' # Define some example edges
#' edges <- data.frame(regulator = c("A", "B", "C"),
#'                     target = c("D", "E", "F"),
#'                     regulatorType = c("miRNA", "TF", "TF"))
#' 
#' # Create TRNet object
#' trn <- TRNet(edges, FALSE, "RNA-seq")
#' 
#' # Generate igraph object
#' generatePlot(trn)
#' 
methods::setGeneric("generatePlot",
                    function(trn) {
                      standardGeneric("generatePlot")
                    })
methods::setMethod("generatePlot", "TRNet",
  function(trn) {
    # Parse vertex metadata
    vertexMetadata <- parseVertexMetadata(trn)
    # Generate igraph object
    network <- igraph::graph_from_data_frame(trn@TRNmetadata,
                                             directed = TRUE,
                                             vertices = vertexMetadata)
  })


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
#' @export
#' 
#' @examples 
#' # Define some example edges
#' edges <- data.frame(regulator = c("A", "B", "C"),
#'                    target = c("D", "E", "F"),
#'                    regulatorType = c("miRNA", "TF", "TF"))
#' 
#' # Create TRNet object
#' trn <- TRNet(edges, FALSE, "RNA-seq")
#' 
#' # Parse vertex metadata
#' parseVertexMetadata(trn)
#' 
methods::setGeneric("parseVertexMetadata",
                    function(trn) {
                      standardGeneric("parseVertexMetadata")
                    })
methods::setMethod("parseVertexMetadata", "TRNet",
  function(trn) {
    edgeMetadata <- trn@TRNmetadata
    # Target vertices
    targetVertices <- data.frame(name = unique(vertexMetadata$target),
                                type = "gene")
    # Regulator vertices
    regulatorVertices <- edgeMetadata[, c("regulator", "regulatorType")]
    colnames(regulatorVertices) <- c("name", "type")
    regulatorVertices <- regulatorVertices[!duplicated(regulatorVertices), ]
    # Combine target and regulator vertices
    vertexMetadata <- rbind(targetVertices, regulatorVertices)
    return(vertexMetadata)
  })