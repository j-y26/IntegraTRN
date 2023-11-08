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
#' 