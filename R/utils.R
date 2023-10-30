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
