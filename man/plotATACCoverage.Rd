% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotDEOmics.R
\name{plotATACCoverage}
\alias{plotATACCoverage}
\title{Peak coverage plot for differential accessible regions}
\usage{
plotATACCoverage(objMOList, title = "ATAC Peaks over Chromosomes")
}
\arguments{
\item{objMOList}{An MOList object containing the differential accessible
regions. The objMOList object must contain a DEATAC
element, which is a result of providing the ATACseq data
and performing the diffOmics() function.}

\item{title}{The title for the plot, default is "ATAC Peaks over Chromosomes"}
}
\value{
A ggplot object
}
\description{
This function generates a peak coverage plot for differential
             accessible regions based on the ATACseq data. This is a wrapper
             function for the \code{\link[ChIPseeker]{covplot}} function in
             the ChIPseeker package.
}
\examples{
\dontrun{
# Use the package-provided example data
data("expMOList")

# Annotates the peaks in the example data and performs motif enrichment
# Using human annotations as an example
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")
library("BSgenome.Hsapiens.UCSC.hg38")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annoDb <- "org.Hs.eg.db"
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

# Load PWMs from JASPAR
data("jasparVertebratePWM")

# Performs the annotation and motif enrichment
expMOList <- annotateATACPeaksMotif(expMOList,
  tssRegion = c(-3000, 3000),
  TxDb = txdb,
  annoDb = annoDb,
  bsgenome = bsgenome,
  pwmL = jasparVertebratePWM
)

# Plotting the coverage plot
plotATACCoverage(expMOList)
}

}
\references{
\insertRef{yu2015chipseeker}{IntegraTRN}

\insertRef{lawrence2013software}{IntegraTRN}
}
