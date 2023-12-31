% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TRNet.R
\name{TRNet}
\alias{TRNet}
\title{TRNet constructor}
\usage{
TRNet(TRNmetadata, predicted, omics)
}
\arguments{
\item{TRNmetadata}{A data frame containing all interaction data for the
transcriptional regulatory network in adjacency list
format. The columns minimally includes "regulator",
"target", and "regulatorType".}

\item{predicted}{A logical value indicating whether the transcriptional
regulatory network contains inferred interactions.}

\item{omics}{A character string indicating the omics data used to construct
the transcriptional regulatory network.}
}
\value{
A TRNet object
}
\description{
This function is a constructor of the TRNet object
}
\examples{
# Define some example edges
edges <- data.frame(
  regulator = c("A", "B", "C"),
  target = c("D", "E", "F"),
  regulatorType = c("miRNA", "TF", "TF")
)

# Create TRNet object
trn <- TRNet(edges, FALSE, "RNA-seq")

}
\references{
\insertRef{csardi2006igraph}{IntegraTRN}
}
