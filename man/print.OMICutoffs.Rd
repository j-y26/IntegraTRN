% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructNetwork.R
\name{print.OMICutoffs}
\alias{print.OMICutoffs}
\title{Print method for the OMICutoffs object}
\usage{
\method{print}{OMICutoffs}(x, ...)
}
\arguments{
\item{x}{An OMICutoffs object}

\item{...}{Other arguments passed to the print function}
}
\value{
The cutoffs for the omics data
}
\description{
This function prints the cutoffs for the omics data to the
             console.
}
\examples{
omiCutoffs <- setOmicCutoffs()
print(omiCutoffs)

}
