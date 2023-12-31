% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MOList.R
\name{newMOList}
\alias{newMOList}
\title{Construct a new MOList object}
\usage{
newMOList(
  RNAseq,
  RNAGroupBy,
  smallRNAseq,
  smallRNAGroupBy,
  proteomics,
  proteomicsGroupBy,
  peakCond1,
  peakCond2
)
}
\arguments{
\item{RNAseq}{A numeric matrix containing the RNAseq data}

\item{RNAGroupBy}{A vector of grouping information for the RNAseq data. Must
be a vector of the same length as the number of samples in
the RNAseq data}

\item{smallRNAseq}{A numeric matrix containing the smallRNAseq data}

\item{smallRNAGroupBy}{A vector of grouping information for the smallRNAseq
data. Must be a vector of the same length as the
number of samples in the smallRNAseq data}

\item{proteomics}{A numeric matrix containing the proteomics data}

\item{proteomicsGroupBy}{A vector of grouping information for the proteomics
data. Must be a vector of the same length as the
number of samples in the proteomics data}

\item{peakCond1}{A data frame containing the differentially accessible ATAC
peaks for condition 1. Should follow the BED format}

\item{peakCond2}{A data frame containing the differentially accessible ATAC
peaks for condition 2. Should follow the BED format}
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
This function constructs a new MOList object, which is used to
             store the multi-omics data. Fulfills the CREATION mode of the
             constructor function.
}
\keyword{internal}
