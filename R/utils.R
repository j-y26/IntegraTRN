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