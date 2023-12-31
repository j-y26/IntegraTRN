% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diffOmicsAnalysis.R
\name{diffExprEdgeR}
\alias{diffExprEdgeR}
\title{Perform differential expression analysis on the count-based omics data using
EdgeR}
\usage{
diffExprEdgeR(filteredCounts, groupBy, batch = NULL)
}
\arguments{
\item{filteredCounts}{A numeric matrix containing the filtered count-based
omics data}

\item{groupBy}{A vector specifying the grouping information for the omics
data, must be the same length as the number of samples in the
omics data}

\item{batch}{A character vector specifying the batch information for the
omics data, must be the same length as the number of samples in
the omics data, used for batch correction. Can be NULL if no
batch correction is needed.}
}
\value{
A DETag object containing the differential expression analysis
        results, and the method DESeq2
\itemize{
\item \code{DEResult}: A data frame containing the differential expression
                       analysis results
\item \code{method}: The character string "EdgeR"
}
}
\description{
This function performs differential expression analysis on the
             count-based omics data using EdgeR. The RNAseq, small RNAseq,
             and protein data are supported for differential expression
             analysis.
}
\references{
\insertRef{robinson2010edger}{IntegraTRN}
}
\keyword{internal}
