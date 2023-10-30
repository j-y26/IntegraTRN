# Purpose: Performing differential analysis on the omics data
# Author: Jielin Yang
# Date: 2023-10-29
# Version: 1.0
# Bugs and Issues: None

# Validate the input data and annotations on the samples
# @param objMOList A MOList object containing the omics data
# @param annoList A list containing the annotations on the samples of the omics
#                 data, must in the order of RNAseq, small RNAseq, and protein
validateDataAnno <- function(objMOList, annoList) {
  # Validate the correctness of the MOList object
  validateMOList(objMOList)
  # Check of annotations on the samples
  if (matchVecToMatrix(annoList$rnaseq, objMOList@RNAseq) ||
    matchVecToMatrix(annoList$smallRna, objMOList@smallRNAseq) ||
    matchVecToMatrix(annoList$protein, objMOList@protein)) {
    stop("The annotations on the samples do not match the omics data")
  } else {
    # Do nothing
  }
}

#' Perform differential analysis on the omics data
#'
#' @description This function performs data filtering, normalization, batch
#'              correction, and differential analysis on each of the omics data.
#'              RNAseq, small RNAseq, and protein data are supported for
#'              differential expression analysis. The ATACseq data are used
#'              used to calculate the differential accessible regions.
#' @param objMOList A MOList object containing the omics data
#' @param rnaseqBatch A character vector specifying the batch information for
#'                    the RNA-seq data, must be the same length as the number
#'                    of samples in the RNA-seq data, used for batch correction
#' @param smallRnaBatch A character vector specifying the batch information for
#'                      the small RNA data, must be the same length as the
#'                      number of samples in the small RNAseq data, used for
#'                      batch correction
#' @param proteinBatch A character vector specifying the batch information for
#'                     the protein data, must be the same length as the number
#'                     of samples in the protein data, used for batch correction
#' @return An MOList object containing the differential analysis results
#' @export
diffOmicsAnalysis <- function(objMOList,
                              rnaseqBatch = NULL,
                              smallRnaBatch = NULL,
                              proteinBatch = NULL) {
  # Validating inputs
  validateDataAnno(objMOList, list(
    rnaseq = rnaseqBatch,
    smallRna = smallRnaBatch,
    protein = proteinBatch
  ))
}
