# Purpose: Performing differential analysis on the omics data
# Author: Jielin Yang
# Date: 2023-10-29
# Version: 1.0
# Bugs and Issues: None


# Define a global variable for the count-based omics data
COUNT_OMICS <- c("RNAseq", "smallRNAseq", "proteomics")


#' Validate the input data and annotations on the samples
#'
#' @keywords internal
#'
#' @description This function validates the input data and annotations on the
#'             samples. The input data must be a MOList object, and the
#'            annotations on the samples must be a list containing the
#'           annotations on the RNAseq, small RNAseq, and protein data, in the
#'         order of RNAseq, small RNAseq, and protein.
#'
#' @param objMOList A MOList object containing the omics data
#' @param annoList A list containing the annotations on the samples of the omics
#'                 data, must contain fields RNAseq, smallRNAseq, and proteomics
#' \itemize{
#'  \item RNAseq: A character/numeric vector for annotating each sample in the
#'                RNAseq data, must be the same length as the number of samples
#'                in the RNAseq data
#' \item smallRNAseq: A character/numeric vector for annotating each sample in
#'                    the small RNAseq data, must be the same length as the
#'                    number of samples in the small RNAseq data
#' \item proteomics: A character/numeric vector for annotating each sample in
#'                   the protein data, must be the same length as the number of
#'                   samples in the protein data
#' }
#'
#' @return NULL
#'
#' @examples
#' \dontrun{
#' # Assuming myMOList is a MOList object
#'
#' # Create an example annotation list for batch correction
#' annoList <- list(
#'   RNAseq = rep(c("A", "B"), each = 5),
#'   smallRNAseq = rep(c("A", "B"), each = 3),
#'   proteomics = rep(c("A", "B"), each = 3)
#' )
#'
#' # Validate the input data and annotations
#' validateDataAnno(myMOList, annoList)
#' }
#'
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
  return(invisible(NULL))
}


#' Filter the gene counts for the count-based omics data
#'
#' @description Filtering is based on the design of the experiment. If the
#'              samples are only grouped into 2 conditions, then the genes with
#'              counts in more than 1 counts per million (CPM) in at least
#'              min(#samples in condition 1, #samples in condition 2) samples
#'              are kept. If the samples are grouped by a continuous variable,
#'              then the genes with counts in more than 1 CPM in at least
#'              30% of the samples are kept.
#'
#' @param objMOList A MOList object containing the omics data
#' @param omic A character string specifying the omics data to be filtered
#'        must be one of "RNAseq", "smallRNAseq", and "proteomics"
#'
#' @return An MOList object containing the filtered omics data
#'
#' @export
#' @importFrom edgeR cpm
#'
#' @references
#' Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for
#' differential expression analysis of digital gene expression data.
#' Bioinformatics. 2010 Jan 1;26(1):139-40. doi: 10.1093/bioinformatics/btp616.
#' Epub 2009 Nov 11. PMID: 19910308; PMCID: PMC2796818.
#'
#' @examples
#' # Create example RNAseq data
#' rnaseq <- matrix(sample(0:100, 1000, replace = TRUE), nrow = 100, ncol = 10)
#' rnaseq[1:20, 1:8] <- 0
#' group <- rep(c("A", "B"), each = 5)
#'
#' # Create example MOList object
#' objMOList <- MOList(RNAseq = rnaseq, RNAGroupBy = group)
#'
#' # Before filtering, check the dimensions of the RNAseq data
#' dim(getCounts(objMOList, "RNAseq"))
#'
#' # Filter the RNAseq data
#' objMOList <- filterGeneCounts(objMOList, "RNAseq")
#'
#' # After filtering, check the dimensions of the RNAseq data
#' dim(getCounts(objMOList, "RNAseq")) # should be lower than the original
#'
filterGeneCounts <- function(objMOList, omic) {
  if (!omics %in% COUNT_OMICS) {
    stop("The input omics data is not supported for filtering")
  }

  omicData <- getCounts(objMOList, omic)
  if (is.null(omicData)) {
    # Nothing to do, return the original object
    return(objMOList)
  } else {
    # Filter the omics data
    # Convert the omics data into counts per million (CPM)
    cpmData <- edgeR::cpm(omicData)

    # Use the grouping information to determine the filtering threshold
    groupBy <- switch(omic,
      RNAseq = objMOList@RNAseqSamples$groupBy,
      smallRNAseq = objMOList@smallRNAseqSamples$groupBy,
      proteomics = objMOList@proteomicsSamples$groupBy
    )
    if (length(unique(groupBy)) == 2) { # Two level
      minSamples <- min(table(groupBy))
      sel <- rowSums(cpmData > 1) >= minSamples
      omicData <- omicData[sel, ]
    } else { # Continuous variable
      sel <- rowSums(cpmData > 1) >= 0.3 * ncol(omicData)
      omicData <- omicData[sel, ]
    }

    # Free up memory for cpm, because it is a large matrix
    rm(cpmData)

    # Update the MOList object
    switch(omic,
      RNAseq = RNAseq(objMOList) <- omicData,
      smallRNAseq = smallRNAseq(objMOList) <- omicData,
      proteomics = protein(objMOList) <- omicData
    )
    return(objMOList)
  }
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
diffOmics <- function(objMOList,
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
  for (omics in COUNT_OMICS) {
    objMOList <- filterGeneCounts(objMOList, omics)
  }
}
