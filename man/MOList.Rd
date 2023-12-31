% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MOList.R
\name{MOList}
\alias{MOList}
\title{Constructor for the MOList object}
\usage{
MOList(
  objMOList = NULL,
  RNAseq = NULL,
  RNAGroupBy = NULL,
  smallRNAseq = NULL,
  smallRNAGroupBy = NULL,
  proteomics = NULL,
  proteomicsGroupBy = NULL,
  ATACpeak1 = NULL,
  ATACpeak2 = NULL
)
}
\arguments{
\item{objMOList}{An object of class MOList for appending/exchanging omics
data. If provided, the MOList object will be modified by
appending/exchanging the omics data. If not provided, a new
MOList object will be constructed. Provision of the
objMOList parameter determines the mode of the constructor
function}

\item{RNAseq}{A numeric matrix containing the RNAseq data, with row names as
gene names and column names as sample names}

\item{RNAGroupBy}{A vector of grouping information for the RNAseq data, used
to perform differential expression analysis. Must be a
vector of the same length as the number of samples in the
RNAseq data}

\item{smallRNAseq}{A numeric matrix containing the smallRNAseq data, with row
names as gene names and column names as sample names}

\item{smallRNAGroupBy}{A vector of grouping information for the smallRNAseq
data, used to perform differential expression analysis.
 Must be a vector of the same length as the number of
 samples in the smallRNAseq data}

\item{proteomics}{A numeric matrix containing the proteomics data, with row
names as gene names and column names as sample names}

\item{proteomicsGroupBy}{A vector of grouping information for the proteomics
data, used to perform differential expression
analysis. Must be a vector of the same length as
the number of samples in the proteomics data}

\item{ATACpeak1}{A character string containing the path to the BED file
containing the unified ATAC peaks for condition 1. The
BED file should NOT contain a header line. Alternatively,
provide a data frame containing the ATAC peaks for condition
with format consistent with the BED definition}

\item{ATACpeak2}{A character string containing the path to the BED file
containing the unified ATAC peaks for condition 2. The
BED file should NOT contain a header line. Alternatively,
provide a data frame containing the ATAC peaks for condition
with format consistent with the BED definition}
}
\value{
An object of class MOList
}
\description{
This function is a constructor of the MOList object, which is
             used to store the multi-omics data. The MOList object can be
             constructed by providing the RNAseq, smallRNAseq,
             proteomics, and ATAC peaks data at once, or appending/exchanging
             each omics data separately. The RNAseq data is required for the
             initial construction of the MOList object. If data of existing
             omics type is provided, old data will be replaced.
}
\note{
The users should ensure that each BED file used as input
                contains the chromosome regions that are found to have
                increased accessibility in each condition. The BED regions
                can be constructed either by merging the peaks from all
                replicates in each condition using a variable-length method
                or a fixed-width approach, as indicated in Grandi et al. 2022
}
\examples{
# Generating some example data
# Note that the exact value could differ based on the random seed
RNAseq <- matrix(sample(1:100, 100, replace = TRUE), ncol = 10)
colnames(RNAseq) <- paste0("sample_", seq_len(ncol(RNAseq)))
rownames(RNAseq) <- paste0("gene_", seq_len(nrow(RNAseq)))
RNAGroupBy <- rep(c("A", "B"), each = 5)

smallRNAseq <- matrix(sample(1:100, 20, replace = TRUE), ncol = 4)
smallRNAGroupBy <- rep(c("A", "B"), each = 2)

proteomics <- matrix(sample(1:100, 30, replace = TRUE), ncol = 6)
proteomicsGroupBy <- rep(c("A", "B"), each = 3)

# Example 1: Constructing the MOList object at once

# Constructing the MOList object at once
objMOList1 <- MOList(
  RNAseq = RNAseq,
  RNAGroupBy = RNAGroupBy,
  smallRNAseq = smallRNAseq,
  smallRNAGroupBy = smallRNAGroupBy,
  proteomics = proteomics,
  proteomicsGroupBy = proteomicsGroupBy
)

# Example 2: Constructing the MOList object step by step

# The MOList object can be constructed minimally with the RNAseq data
# Further omics data can be appended/exchanged later
objMOList2 <- MOList(RNAseq = RNAseq, RNAGroupBy = RNAGroupBy)
RNAseq2 <- matrix(sample(1:100, 100, replace = TRUE), ncol = 10)
colnames(RNAseq2) <- paste0("sample_", seq_len(ncol(RNAseq2)))
objMOList2 <- MOList(objMOList2, RNAseq = RNAseq2, RNAGroupBy = RNAGroupBy)
objMOList2 <- MOList(objMOList2,
  proteomics = proteomics,
  proteomicsGroupBy = proteomicsGroupBy
)

# Example 3: Exchanging some inputs for the MOList object

# The MOList object can be modified by exchanging some inputs
objMOList3 <- objMOList2
RNAseq3 <- RNAseq
RNAseq3[1:10, 1:5] <- 0
objMOList3 <- MOList(objMOList3, RNAseq = RNAseq3, RNAGroupBy = RNAGroupBy)

}
\references{
\insertRef{genomictools}{IntegraTRN}
}
