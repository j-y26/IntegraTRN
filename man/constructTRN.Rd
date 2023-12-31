% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructNetwork.R
\name{constructTRN}
\alias{constructTRN}
\title{Construct a transcriptional regulatory network}
\usage{
constructTRN(
  objMOList,
  omiCutoffs,
  useOmics = c("RNAseq", "smallRNAseq", "proteomics", "ATACseq", "extInteractions"),
  smallRNAtypes = "all",
  targetDirection = c("up", "down", "both"),
  predicted = TRUE,
  ntree = 1000,
  nthreads = 1,
  treeMethod = "RF",
  seed = 91711
)
}
\arguments{
\item{objMOList}{A MOList object containing all omics data that the user
wants to use to construct the network}

\item{omiCutoffs}{A OMICutoffs object containing the cutoffs for the omics
data}

\item{useOmics}{A character vector containing the omics data that the user
wants to use to construct the network. The available omics
data are "RNAseq", "smallRNAseq", "proteomics", "ATACseq",
and "extInteractions". By default, all the omics data will
be used.}

\item{smallRNAtypes}{A character vector containing the small RNA types that
the user wants to use to construct the network. The
available types are "miRNA", "piRNA", "snRNA", "snoRNA",
"tRNA", and "circRNA". If "all" is specified, then all
the small RNA types will be used.}

\item{targetDirection}{A single-length character vector indicating the
direction of expression change of the target genes.
Default to "both". Other options include "up" and
"down".}

\item{predicted}{A logical value indicating whether to perform predicted
inference of small RNA - mRNA interactions. Default to TRUE.
Prediction utilizes the correlation of expression between
small RNAs and mRNAs, which is highly recommended for adding
an additional level of condition-specificity to the network.}

\item{ntree}{The number of trees to grow in the random forest model, Must be
a positive integer. Default is 1000.}

\item{nthreads}{The number of threads to use for parallel processing. Must be
a positive integer between 1 and the number of available
cores on the user's machine. Default is 1.}

\item{treeMethod}{The method to use for the random forest model. Either "RF"
or "ET" for random forest and extra trees, respectively.
See the documentation for the \code{\link[GENIE3]{GENIE3}}
package for details. Default is "RF".}

\item{seed}{The seed to use for the random forest model. Must be a positive
integer. Default is 91711.}
}
\value{
A TRNet object containing the transcriptional regulatory network
}
\description{
This function constructs a transcriptional regulatory network
             based on the omics data and the external interaction data
             provided. The network is stored as an igraph object in the
             TRNet S4 class, which is a mRNA-TF-smallRNA network.
}
\examples{
# Use the package data
data("expMOList")

# Set the cutoffs for the omics data
omiCutoffs <- setOmicCutoffs() # use default cutoffs

# Construct the network
\dontrun{
myTRNet <- constructTRN(expMOList, omiCutoffs, targetDirection = "up")
plotNetwork(myTRNet)
}

}
\references{
\insertRef{dplyr}{IntegraTRN}

\insertRef{csardi2006igraph}{IntegraTRN}

\insertRef{huynh2010inferring}{IntegraTRN}
}
