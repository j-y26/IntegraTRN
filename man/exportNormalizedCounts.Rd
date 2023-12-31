% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DETag.R
\name{exportNormalizedCounts}
\alias{exportNormalizedCounts}
\alias{exportNormalizedCounts,DETag-method}
\title{Export the normalized counts from a DETag object}
\usage{
exportNormalizedCounts(x)
}
\arguments{
\item{x}{An object of the DETag class}
}
\value{
A matrix containing the normalized counts for the genes, with each
        row representing a gene and each column representing a sample
}
\description{
This function extracts the normalized counts from the DETag
             object, but this only applies to the count-based differential
             analysis
}
\examples{
# Use the package-provided example data
data("expMOList")
deTag <- expMOList$DERNAseq

# Export the normalized counts
exportNormalizedCounts(deTag)

}
