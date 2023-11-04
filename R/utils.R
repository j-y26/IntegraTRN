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


#' Defining DESeq2 colData and design from vectors
#'
#' @keywords internal
#'
#' @param groupBy A vector of group names
#' @param batch A vector of batch names
#'
#' @return A list containing the colData and design, where colData is a data
#'         frame and design is a formula
#'
DESeqDesign <- function(groupBy, batch = NULL) {
  colData <- data.frame(group = groupBy)
  if (!is.null(batch)) {
    if (length(batch) != length(groupBy)) {
      stop("Length of batch and groupBy vectors do not match")
    } else {
      colData$batch <- batch
      design <- ~ batch + group
    }
  } else {
    design <- ~group
  }
  return(list(colData = colData, design = design))
}

#' Transform the result of MatchIt::matchit() to a data frame
#'
#' @keywords internal
#'
#' @param matchResult The result of MatchIt::matchit()
#'
#' @return A data frame containing the matched pairs
#'
matchResultToDF <- function(matchResult) {
  if (is.null(matchResult)) {
    stop("Error generating matches.")
  }
  matchResult <- matchResult$match.matrix
  sampleMatch <- data.frame()
  if (all(grepl(SRNA_SUFFIX, rownames(matchResult)))) {
    sampleMatch$smallRNA <- gsub(SRNA_SUFFIX, "", rownames(matchResult))
    sampleMatch$RNA <- matchResult[, 1]
  } else {
    sampleMatch$smallRNA <- gsub(SRNA_SUFFIX, "", matchResult[, 1])
    sampleMatch$RNA <- rownames(matchResult)
  }
  return(sampleMatch)
}


#' Label the samples with smaller sample size as 1
#'
#' @keywords internal
#'
#' @param sampleDF A data frame containing the sample names as row names
#' @param identifier A string to identify the samples in sample names
#' @param colname The column name to be added to the data frame
#'
#' @return A data frame with the column added
#'
labelSmallSizeGroup <- function(sampleDF, identifier, colname) {
  identified <- grepl(identifier, row.names(sampleDF))
  nonIdentified <- !identified
  if (sum(identified) < sum(nonIdentified)) {
    sampleDF[identified, colname] <- 1
    sampleDF[nonIdentified, colname] <- 0
  } else {
    sampleDF[identified, colname] <- 0
    sampleDF[nonIdentified, colname] <- 1
  }
  return(sampleDF)
}


# [END]
