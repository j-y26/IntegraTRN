% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MOList.R
\name{getSampleInfo}
\alias{getSampleInfo}
\title{Retrieving the sample information from the MOList object}
\usage{
getSampleInfo(x, experiment)
}
\arguments{
\item{x}{A MOList object containing the omics data}

\item{experiment}{A character string specifying the experiment type, must
be one of "RNAseq", "smallRNAseq", and "proteomics"}
}
\value{
A list containing the sample names and grouping information for the
        specified omics data
}
\description{
Retrieving the sample information from the MOList object
}
\keyword{internal}
