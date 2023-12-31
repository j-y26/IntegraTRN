% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MOList.R
\name{RNAseq<-}
\alias{RNAseq<-}
\title{Setter for the RNAseq slot}
\usage{
RNAseq(x) <- value
}
\arguments{
\item{x}{An object of class MOList for appending/exchanging omics data}

\item{value}{A numeric matrix containing the RNAseq data}
}
\value{
An object of class MOList
\itemize{
\item \code{RNAseq}: A numeric matrix containing the RNAseq data
\item \code{RNAseqSamples}: A list containing the sample names and grouping
                        information for the RNAseq data
\item \code{smallRNAseq}: A numeric matrix containing the smallRNAseq data
\item \code{smallRNAseqSamples}: A list containing the sample names and
                             grouping information for the smallRNAseq data
\item \code{proteomics}: A numeric matrix containing the proteomics data
\item \code{proteomicsSamples}: A list containing the sample names and
                          grouping information for the proteomics data
\item \code{ATACpeaks}: A list containing the differentially accessible ATAC
                        peaks for condition 1 and condition 2
}
}
\description{
This function is a setter for the RNAseq slot of the MOList.
             The user are NOT allowed to use this function directly, since
             use without caution may cause the MOList object to be invalid.
}
\keyword{internal}
