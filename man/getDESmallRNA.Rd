% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processSmallRNA.R
\name{getDESmallRNA}
\alias{getDESmallRNA}
\title{Retrieve differentially expressed small RNAs with specified types}
\usage{
getDESmallRNA(
  objMOList,
  padj = 0.05,
  log2fc = 0,
  type = c("miRNA", "piRNA", "tRNA", "circRNA", "snRNA", "snoRNA")
)
}
\arguments{
\item{objMOList}{A object of class MOList}

\item{padj}{A numeric value indicating the adjusted p-value threshold.
Default is 0.05.}

\item{log2fc}{A numeric value indicating the log2 fold change threshold.
Default is 0.}

\item{type}{A character vector indicating the small RNA types to be
retrieved. The default is all small RNA types.}
}
\value{
A data frame containing the differentially expressed small RNAs
        that meet the criteria.
}
\description{
This function retrieves differentially expressed small RNAs
             based on user-specified criteria. The criteria include the
             following: 1. adjusted p-value; 2. log2 fold change; 3. small
             RNA type. The function will return a data frame containing the
             differentially expressed small RNAs that meet the criteria.
}
\examples{
# Example 1: Retrieve all differentially expressed small RNAs
# Use the package provided example data
data(expMOList)

# Retrieve all differentially expressed small RNAs
getDESmallRNA(expMOList)
# or alternatively, specify all the parameters
getDESmallRNA(expMOList,
  padj = 0.05,
  log2fc = 0,
  type = c(
    "miRNA", "piRNA", "tRNA",
    "circRNA", "snRNA", "snoRNA"
  )
)

# Example 2: Retrieve differentially expressed miRNAs with default cutoffs
getDESmallRNA(expMOList, type = "miRNA")

}
\references{
\insertRef{dplyr}{IntegraTRN}
}
