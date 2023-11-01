# Purpose: Utility functions more streamlined code
# Author: Jielin Yang
# Date: 2023-10-29
# Version: 1.0
# Bugs and Issues: None


#' Identify the sample names from a count matrix
#'
#' @keywords internal
#'
#' @param countMatrix A count matrix with samples as columns and features as
#'                    rows
#'
#' @return A vector of sample names
#'
#' @examples
#' # Create example count matrix
#' countMatrix <- matrix(sample(0:100, 1000, replace = TRUE),
#'   nrow = 100, ncol = 10
#' )
#'
#' # Get sample names
#' getSampleNames(countMatrix)
#'
getSampleNames <- function(countMatrix) {
  if (is.null(colnames(countMatrix))) {
    sampleNames <- paste0("sample_", seq_len(ncol(countMatrix)))
  } else {
    sampleNames <- colnames(countMatrix)
  }
  return(sampleNames)
}


#' Match vector length to matrix column length
#'
#' @keywords internal
#'
#' @param vector A vector
#' @param matrix A matrix
#'
#' @return A boolean indicating whether the vector length matches the matrix
#'         column length
#' @examples
#' # Create example vector and matrix
#' vec <- 1:10
#' mat <- matrix(1:100, nrow = 10, ncol = 10)
#'
#' # Match vector length to matrix column length
#' matchVecToMatrix(vec, mat)
#'
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


# [END]
