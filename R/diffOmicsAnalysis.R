# Purpose: Performing differential analysis on the omics data
# Author: Jielin Yang
# Date: 2023-10-29
# Version: 1.0
# Bugs and Issues: None


# Validate the input data and annotations on the samples
# Level: private
# @param objMOList A MOList object containing the omics data
# @param annoList A list containing the annotations on the samples of the omics
#                 data, must in the order of RNAseq, small RNAseq, and protein
validateDataAnno <- function(objMOList, annoList) {
  # Validate the correctness of the MOList object
  validateMOList(objMOList)
  # Check of annotations on the samples
  if (matchVecToMatrix(annoList$RNAseq, objMOList@RNAseq) ||
    matchVecToMatrix(annoList$smallRNAseq, objMOList@smallRNAseq) ||
    matchVecToMatrix(annoList$proteomics, objMOList@protein)) {
    stop("The annotations on the samples do not match the omics data")
  } else {
    # Do nothing
  }
}


# Filter the gene counts for the count-based omics data
# Level: private
#
# @description Filtering is based on the design of the experiment. If the
#              samples are only grouped into 2 conditions, then the genes with
#              counts in more than 1 counts per million (CPM) in at least 
#              min(#samples in condition 1, #samples in condition 2) samples
#              are kept. If the samples are grouped by a continuous variable,
#              then the genes with counts in more than 1 CPM in at least
#              30% of the samples are kept.
# @param objMOList A MOList object containing the omics data
# @param omics A character string specifying the omics data to be filtered
#        must be one of "RNAseq", "smallRNAseq", and "proteomics"
# @return An MOList object containing the filtered omics data
filterGeneCounts <- function(objMOList, omic) {
}






#' Perform differential analysis on the omics data
#'
#' @description This function performs data filtering, normalization, batch
#'              correction, and differential analysis on each of the omics data.
#'              RNAseq, small RNAseq, and protein data are supported for
#'              differential expression analysis. The ATACseq data are used
#'              used to calculate the differential accessible regions.
#'              While performing the differential analysis, the original count
#'              data for the RNAseq, small RNAseq, and protein data are
#'              filtered and subsequently normalized by library size. After
#'              calling diffOmicsAnalysis, the users can use the S4 methods
#'              of the MOList class to extract filtered and normalized read
#'              counts for each of the omics data.
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
#' 
diffOmicsAnalysis <- function(objMOList,
                              rnaseqBatch = NULL,
                              smallRnaBatch = NULL,
                              proteinBatch = NULL) {
  # Validating inputs
  validateDataAnno(objMOList, list(
    RNAseq = rnaseqBatch,
    smallRNAseq = smallRnaBatch,
    proteomics = proteinBatch
  ))

  # Data filtering for count-based omics data
  for (omics in c("RNAseq", "smallRNAseq", "proteomics")) {
    objMOList <- filterMOList(objMOList, omics)
  }

}


