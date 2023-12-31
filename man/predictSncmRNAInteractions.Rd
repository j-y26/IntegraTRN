% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructNetwork.R
\name{predictSncmRNAInteractions}
\alias{predictSncmRNAInteractions}
\title{Predict small RNA - mRNA interactions}
\usage{
predictSncmRNAInteractions(
  objMOList,
  omiCutoffs,
  rnaTopTag,
  smallRNATypes = SMALLRNA_CATEGORIES,
  ntree = 1000,
  nthreads = 1,
  treeMethod = "RF",
  seed = 91711
)
}
\arguments{
\item{objMOList}{A MOList object containing the omics data}

\item{omiCutoffs}{A OMICutoffs object containing the cutoffs for the omics
data}

\item{rnaTopTag}{A TOPTag object containing the top differential genes from
the RNAseq data}

\item{smallRNATypes}{A character vector containing the small RNA types that
the user wants to use to predict the interactions. The
available types are "miRNA", "piRNA", "snRNA", "snoRNA",
"tRNA", and "circRNA".}

\item{ntree}{The number of trees to grow in the random forest model}

\item{nthreads}{The number of threads to use for parallel processing}

\item{treeMethod}{The method to use for the random forest model. See the
documentation for the randomForest package for details.}

\item{seed}{The seed to use for the random forest model}
}
\value{
a list of predicted interactions combined for all specified small
        RNA types as the regulators
}
\description{
This function predicts the small RNA - mRNA interactions based
             on the smallRNAseq data and the RNAseq data.
}
\keyword{internal}
