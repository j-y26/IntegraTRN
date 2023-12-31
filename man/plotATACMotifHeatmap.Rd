% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotDEOmics.R
\name{plotATACMotifHeatmap}
\alias{plotATACMotifHeatmap}
\title{Motif heatmap of differential accessible regions of enriched motifs}
\usage{
plotATACMotifHeatmap(
  objMOList,
  pValueAdj = 0.05,
  pValue = NULL,
  log2FEnrich = NULL
)
}
\arguments{
\item{objMOList}{An MOList object containing the differential accessible
regions}

\item{pValueAdj}{The cutoff for adjusted p-value, default is 0.05}

\item{pValue}{The cutoff for p-value. If pValue is specified, pValueAdj
will be ignored}

\item{log2FEnrich}{The cutoff for log2 fold enrichment, default is NULL}
}
\value{
A ComplexHeatmap object
}
\description{
This function generates a combined profile and heatmap of
             differential accessible regions of enriched motifs based on the
             ATACseq data. This is a wrapper function for the
             plotMotifHeatmaps function in the ChIPseeker package.
}
\details{
The motif enrichment analysis can be exceptionally stringent for
         ATACseq data. The users can specify the cutoff for the adjusted
         p-value or the regular p-value. The users can also specify the
         log2 enrichment fold change cutoff. Please use this function to
         explore the enriched motifs in the differential accessible regions
         and decide the best combination of cutoffs for constructing the
         transcriptional regulatory network. If pValue is specified,
         pValueAdj will be ignored.
}
\examples{
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

# Example 1: Plot the motif heatmap by default parameters
# This example should not show up any motif but instead issue an error
# and returns NULL since no motif are found to be differentially
# enriched at adjusted p value < 0.05

\dontrun{
plotATACMotifHeatmap(expMOList)
}

# Example 2: Plot the motif heatmap with unadjusted p-value cutoff of 0.01
plotATACMotifHeatmap(expMOList, pValue = 0.01)

# Example 3: Plot the motif heatmap with log2 fold enrichment cutoff of 1
plotATACMotifHeatmap(expMOList, pValue = 0.01, log2FEnrich = 1)

}
\references{
\insertRef{machlab2022monalisa}{IntegraTRN}

\insertRef{SummarizedExperiment}{IntegraTRN}

\insertRef{gu2016complex}{IntegraTRN}
}
