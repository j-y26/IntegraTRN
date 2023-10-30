# Purpose: Utility functions more streamlined code
# Author: Jielin Yang
# Date: 2023-10-29
# Version: 1.0
# Bugs and Issues: None


# Identify the sample names from a count matrix
# Level: private
# @param countMatrix A count matrix with samples as columns and features as rows
# @return A vector of sample names
getSampleNames <- function(countMatrix) {
  if (is.null(colnames(countMatrix))) {
    sampleNames <- paste0("sample_", seq_len(ncol(countMatrix)))
  } else {
    sampleNames <- colnames(countMatrix)
  }
  return(sampleNames)
}


# Match vector length to matrix column length
# Level: private
# @param vector A vector
# @param matrix A matrix
# @return A boolean indicating whether the vector length matches the matrix
#         column length
matchVecToMatrix <- function(vec, mat) {
  if (is.null(vec) && is.null(mat)) {
    # Do nothing
  } else if (xor(is.null(vec), is.null(mat))) {
    return(FALSE)
  } else if (length(vec) != ncol(mat)) {
    return(FALSE)
  } else if (any(is.na(vec))) {
    return(FALSE)
  }
  return(TRUE)
}

# Retrieving the grouping information for a specific omics data
# Level: private
# @param objMOList A MOList object containing the omics data
# @param experiment A character string specifying the experiment type, must
#                   be one of "RNAseq", "smallRNAseq", and "proteomics"
# @return A vector of grouping information
getGroupingInfo <- function(objMOList, experiment) {
  switch(experiment,
    RNAseq = objMOList@RNAseq$group,
    smallRNAseq = objMOList@smallRNAseq$group,
    proteomics = objMOList@protein$group
  )
}
