# Purpose: Creating the Multi-Omics List S4 class
# Author: Jielin Yang
# Date: 2023-10-29
# Version: 1.0
# Bugs and Issues: None


# Global variables
CHROMINFO <- c("chrom", "chromStart", "chromEnd")


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
  } else if (!is.numeric(matrix) || !is.matrix(matrix)) {
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
  if (any(is.na(peakCond1[, CHROMINFO])) ||
    any(is.na(peakCond2[, CHROMINFO]))) {
    stop("Missing values in the chromosome information in the ATAC peaks.")
  } else {
    # Do nothing
  }
}


# Construct a new MOList object
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
# @return An object of class MOList
newMOList <- function(RNAseq,
                      RNAGroupBy,
                      smallRNAseq,
                      smallRNAGroupBy,
                      proteomics,
                      proteomicsGroupBy,
                      peakCond1 = NA,
                      peakCond2 = NA) {
  objMOList <- new("MOList",
    RNAseq = RNAseq,
    RNAseqSamples = list(
      samples = colnames(RNAseq),
      groupBy = RNAGroupBy
    ),
    smallRNAseq = smallRNAseq,
    smallRNAseqSamples = list(
      samples = colnames(smallRNAseq),
      groupBy = smallRNAGroupBy
    ),
    proteomics = proteomics,
    proteomicsSamples = list(
      samples = colnames(proteomics),
      groupBy = proteomicsGroupBy
    ),
    ATACpeaks = list(
      peaksCond1 = peakCond1,
      peaksCond2 = peakCond2
    )
  )
  return(objMOList)
}


# Append/exchange omics data for the existing MOList object
# Level: Private
# @param objMOList An object of class MOList for appending/exchanging omics data
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
# @return An object of class MOList
modifyMOList <- function(objMOList,
                         RNAseq,
                         RNAGroupBy,
                         smallRNAseq,
                         smallRNAGroupBy,
                         proteomics,
                         proteomicsGroupBy,
                         peakCond1 = NA,
                         peakCond2 = NA) {
  if (!is.na(RNAseq)) {
    objMOList@RNAseq <- RNAseq
    objMOList@RNAseqSamples$samples <- colnames(RNAseq)
    objMOList@RNAseqSamples$groupBy <- RNAGroupBy
  } else {
    # Do nothing
  }
  if (!is.na(smallRNAseq)) {
    objMOList@smallRNAseq <- smallRNAseq
    objMOList@smallRNAseqSamples$samples <- colnames(smallRNAseq)
    objMOList@smallRNAseqSamples$groupBy <- smallRNAGroupBy
  } else {
    # Do nothing
  }
  if (!is.na(proteomics)) {
    objMOList@proteomics <- proteomics
    objMOList@proteomicsSamples$samples <- colnames(proteomics)
    objMOList@proteomicsSamples$groupBy <- proteomicsGroupBy
  } else {
    # Do nothing
  }
  if (!is.na(peakCond1) && !is.na(peakCond2)) {
    objMOList@ATACpeaks$peaksCond1 <- peakCond1
    objMOList@ATACpeaks$peaksCond2 <- peakCond2
  } else {
    # Do nothing
  }
  return(objMOList)
}














#' Constructor for the MOList object
#'
#' @description This function is a constructor of the MOList object, which is
#'              used to store the multi-omics data. The MOList object can be
#'              constructed by providing the RNAseq, smallRNAseq,
#'              proteomics, and ATAC peaks data at once, or appending/exchanging
#'              each omics data separately. The RNAseq data is required for the
#'              initial construction of the MOList object. If data of existing
#'              omics type is provided, old data will be replaced.
#' @param objMOList An object of class MOList for appending/exchanging omics
#'                  data
#' @param RNAseq A numeric matrix containing the RNAseq data, with row names as
#'               gene names and column names as sample names
#' @param RNAGroupBy A vector of grouping information for the RNAseq data, used
#'                   to perform differential expression analysis
#' @param smallRNAseq A numeric matrix containing the smallRNAseq data, with row
#'                    names as gene names and column names as sample names
#' @param smallRNAGroupBy A vector of grouping information for the smallRNAseq
#'                        data, used to perform differential expression analysis
#' @param proteomics A numeric matrix containing the proteomics data, with row
#'                   names as gene names and column names as sample names
#' @param proteomicsGroupBy A vector of grouping information for the proteomics
#'                          data, used to perform differential expression
#'                          analysis
#' @param pathATACpeak1 The path to a BED file containing the ATAC peaks for
#'                      condition 1
#' @param pathATACpeak2 The path to a BED file containing the ATAC peaks for
#'                      condition 2
#'
#' @return An object of class MOList
#' @export MOList
#'
#'
MOList <- function(objMOList = NA,
                   RNAseq = NA_real_,
                   RNAGroupBy = NA,
                   smallRNAseq = NA_real_,
                   smallRNAGroupBy = NA,
                   proteomics = NA_real_,
                   proteomicsGroupBy = NA,
                   pathATACpeak1 = NA,
                   pathATACpeak2 = NA) {
  # One of objMOList and RNAseq must be given, which dictates whether the
  # constructor creates a new object or appends/exchanges omics data
  if (is.na(objMOList) && is.na(RNAseq)) {
    stop("Please provide either an object of class MOList or the RNAseq data.")
  } else {
    # Do nothing
  }

  # Read the ATAC peaks data if provided
  if (!is.na(pathATACpeak1) && !is.na(pathATACpeak2)) {
    peakCond1 <- GenomicTools.fileHandler::importBed(pathATACpeak1)
    peakCond2 <- GenomicTools.fileHandler::importBed(pathATACpeak2)
  } else if (xor(is.na(pathATACpeak1), is.na(pathATACpeak2))) {
    stop("Please provide both ATAC peaks files.")
  } else {
  }

  # Input validation
  validateMOInputs(
    RNAseq, RNAGroupBy, smallRNAseq, smallRNAGroupBy,
    proteomics, proteomicsGroupBy, peakCond1, peakCond2
  )

  # Construct or modify the MOList object
  if (is.na(objMOList)) {
    # Construct a new MOList object
    newObjMOList <- newMOList(
      RNAseq, RNAGroupBy, smallRNAseq, smallRNAGroupBy,
      proteomics, proteomicsGroupBy, peakCond1, peakCond2
    )
  } else {
    # Modify the existing MOList object
    newObjMOList <- modifyMOList(
      objMOList, RNAseq, RNAGroupBy, smallRNAseq,
      smallRNAGroupBy, proteomics, proteomicsGroupBy,
      peakCond1, peakCond2
    )
  }
  return(newObjMOList)
}




# [END]
