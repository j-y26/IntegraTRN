# Purpose: Creating the Multi-Omics List S4 class
# Author: Jielin Yang
# Date: 2023-10-29
# Version: 1.0
# Bugs and Issues: None

# Define the MOList S4 class
# Slots: RNAseq, RNAseqSamples, smallRNAseq, smallRNAseqSamples, proteomics,
#        proteomicsSamples, ATACpeaks
# Inheritance: list
setClass("MOList",
  # Inheritance
  contains = "list",
  # Slots
  slots = list(
    RNAseq = "numeric",
    RNAseqSamples = "list",
    smallRNAseq = "numeric",
    smallRNAseqSamples = "list",
    proteomics = "numeric",
    proteomicsSamples = "list",
    ATACpeaks = "list"
  ),
  # Prototypes
  prototype = prototype(
    RNAseq = NA_real_,
    RNAseqSamples = list(
      samples = NA,
      groupBy = NA
    ),
    smallRNAseq = NA_real_,
    smallRNAseqSamples = list(
      samples = NA,
      groupBy = NA
    ),
    proteomics = NA_real_,
    proteomicsSamples = list(
      samples = NA,
      groupBy = NA
    ),
    ATACpeaks = list(
      peaksCond1 = NA,
      peaksCond2 = NA
    )
  )
)

# Validate matrix and grouping information
# Level: Private
# @param matrix A numeric matrix containing the sequencing data
# @param groupBy A vector of grouping information for the sequencing data
validateMatrix <- function(matrix, groupBy) {
  if (any(is.na(matrix))) {
    stop("The sequencing data contains NA values. Please check the validity of
    your data.")
  } else if (is.na(groupBy)) {
    stop("Please provide grouping information for the sequencing data.")
  } else if (!is.numeric(matrix)) {
    stop("The sequencing data must be a numeric matrix.")
  } else if (ncol(matrix) != length(groupBy)) {
    stop("The number of samples in the sequencing data must match the length of
    the grouping information.")
  } else {
    # Do nothing
  }
}

# Validate the inputs for the MOList constructor
# Level: Private
# @param RNAseq A numeric matrix containing the RNAseq data
# @param RNAGroupBy A vector of grouping information for the RNAseq data
# @param smallRNAseq A numeric matrix containing the smallRNAseq data
# @param smallRNAGroupBy A vector of grouping information for the smallRNAseq
#                        data
# @param proteomics A numeric matrix containing the proteomics data
# @param proteomicsGroupBy A vector of grouping information for the proteomics
#                          data
# @param peakCond1 A data frame containing the ATAC peaks for condition 1
# @param peakCond2 A data frame containing the ATAC peaks for condition 2
validateMOInputs <- function(RNAseq,
                             RNAGroupBy,
                             smallRNAseq,
                             smallRNAGroupBy,
                             proteomics,
                             proteomicsGroupBy,
                             peakCond1,
                             peakCond2) {
  # Validate the RNAseq data
  validateMatrix(RNAseq, RNAGroupBy)
  # Validate the smallRNAseq data
  if (!is.na(smallRNAseq)) {
    validateMatrix(smallRNAseq, smallRNAGroupBy)
  } else {
    # Do nothing
  }
  # Validate the proteomics data
  if (!is.na(proteomics)) {
    validateMatrix(proteomics, proteomicsGroupBy)
  } else {
    # Do nothing
  }
  # Replicates are required for the sequencing data
  if (ncol(RNAseq) < 2 || ncol(smallRNAseq) < 2 || ncol(proteomics) < 2) {
    stop("At least two replicates are required for each sequencing data.")
  } else {
    # Do nothing
  }
  # Check for missing values in the chromosome information in the ATAC peaks
  chromInfo <- c("chrom", "chromStart", "chromEnd")
  if (any(is.na(peakCond1[, chromInfo])) || 
      any(is.na(peakCond2[, chromInfo]))) {
    stop("Missing values in the chromosome information in the ATAC peaks.")
  } else {
    # Do nothing
  }
}









# [END]
