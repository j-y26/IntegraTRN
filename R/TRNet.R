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