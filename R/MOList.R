# Purpose: Creating the Multi-Omics List S4 class
# Author: Jielin Yang
# Date: 2023-10-29
# Version: 1.0
# Bugs and Issues: None


# Global variables
CHROMINFO <- c("chr", "start", "end")
INTERACTION_FIELDS <- c("ID", "Target")


#' Multi-Omics List (MOList) S4 class
#'
#' @description This class is used to store the multi-omics data, including
#'              RNAseq, smallRNAseq, proteomics, and ATAC peaks data. It
#'              inherits the list class, so that it operates like a list with
#'              subset and replacement methods for non-raw data. The raw data
#'              for each omic type are held within slots and should not be
#'              accessed directly.
#'
#' @slot RNAseq A numeric matrix containing the RNAseq data
#' @slot RNAseqSamples A list with two elements: samples and groupBy. The
#'                     samples element is a character vector containing the
#'                     sample names for the RNAseq data. The groupBy element is
#'                     a character vector containing the grouping information
#'                     for the RNAseq data, used to perform differential
#'                     expression analysis. Both elements have the same length as
#'                     the number of samples in the RNAseq data
#' @slot smallRNAseq A numeric matrix containing the smallRNAseq data
#' @slot smallRNAseqSamples A list with two elements: samples and groupBy.
#'                          Has the same format as the RNAseqSamples slot
#' @slot proteomics A numeric matrix containing the proteomics data
#' @slot proteomicsSamples A list with two elements: samples and groupBy.
#'                         Has the same format as the RNAseqSamples slot
#' @slot ATACpeaks A list containing the ATAC peaks for condition 1 and
#'                 condition 2 that are found to be differentially accessible.
#'                 The users should ensure that each BED file used as input
#'                 contains the chromosome regions that are found to have
#'                 increased accessibility in each condition. The BED regions
#'                 can be constructed either by merging the peaks from all
#'                 replicates in each condition using a variable-length method
#'                 or a fixed-width approach, as indicated in Grandi et al. 2022
#' @slot .Data A list containing differential expression analysis data, used
#'             in the same way as a list
#'
#' @exportClass MOList
#'
#' @references
#' \insertRef{grandi2022chromatin}{IntegraTRN}
#'
setClass("MOList",
  # Inheritance
  contains = "list",
  # Slots
  slots = list(
    RNAseq = "matrix",
    RNAseqSamples = "list",
    smallRNAseq = "matrix",
    smallRNAseqSamples = "list",
    proteomics = "matrix",
    proteomicsSamples = "list",
    ATACpeaks = "list"
  ),
  # Prototypes
  prototype = prototype(
    RNAseqSamples = list(
      samples = NULL,
      groupBy = NULL
    ),
    smallRNAseqSamples = list(
      samples = NULL,
      groupBy = NULL
    ),
    proteomicsSamples = list(
      samples = NULL,
      groupBy = NULL
    ),
    ATACpeaks = list(
      peaksCond1 = NULL,
      peaksCond2 = NULL
    )
  )
)


#' Validate matrix and grouping information
#'
#' @keywords internal
#'
#' @description This function validates the matrix and grouping information for
#'              the sequencing data. The matrix must be a numeric matrix, and
#'              the grouping information must be a vector of the same length as
#'              the number of samples in the sequencing data. Assumes that the
#'              input matrix is not NULL.
#'
#' @param matrix A numeric matrix containing the sequencing data
#' @param groupBy A vector of grouping information for the sequencing data
#'
#' @return NULL
#'
validateMatrix <- function(matrix, groupBy) {
  if (is.null(groupBy)) {
    stop("Grouping information must be provided for provided count data.")
  } else if (any(is.na(matrix))) {
    stop("The sequencing data contains NA values. Please check the validity of
    your data.")
  } else if (any(is.na(groupBy))) {
    stop("Please provide correct grouping information for the count data.")
  } else if (!is.numeric(matrix) || !is.matrix(matrix)) {
    stop("The sequencing data must be a numeric matrix.")
  } else if (!all(matrix >= 0)) {
    stop("The count values in the sequencing data must be non-negative.")
  } else if (!all(matrix == round(matrix))) {
    stop(paste0(
      "The raw input data must be un-normalized counts as obtained ",
      "from the raw counting of the sequencing reads. All count values are ",
      "expected to be integers."
    ))
  } else if (ncol(matrix) != length(groupBy)) {
    stop("The number of samples in the sequencing data must match the length of
    the grouping information.")
  } else if (length(unique(groupBy)) < 2) {
    stop("A minimum of two groups are required for the count data.")
  } else {
    # Do nothing
  }
  return(invisible(NULL))
}


#' Validate BED files for peak information
#'
#' @keywords internal
#'
#' @description This function validates the BED files for peak information.
#'              The BED files are expected to contain minimally the chromosome,
#'              start, and end coordinates of the peaks. This function provides
#'              extendibility for future implementation of additional data
#'              fields in the BED files, such as ChIP-seq peaks. Assumes that
#'              the input peak file is not NULL
#'
#' @param peak A data frame containing the peak information in the BED format
#'
#' @return NULL
#'
#' @references
#' \insertRef{grandi2022chromatin}{IntegraTRN}
#'
validatePeaks <- function(peak) {
  if (!is.data.frame(peak)) {
    stop("The peak information must be able to be read as a data frame.")
  } else if (any(is.na(peak[, CHROMINFO]))) {
    stop("Missing values in the chromosome information in the ATAC peaks.")
  } else if (any(peak[, "start"] < 0) || any(peak[, "end"] < 0)) {
    stop("The start and end coordinates of the peaks must be non-negative.")
  } else if (any(peak[, "start"] > peak[, "end"])) {
    stop("The start coordinates must be smaller than the end coordinates.")
  } else {
    # Do nothing
  }
}


#' Validate the inputs for the MOList constructor
#'
#' @keywords internal
#'
#' @description This function validates the inputs for the MOList constructor.
#'              It uses the validateMatrix function to validate the sequencing
#'              data, and checks for missing values in the chromosome
#'              information in the ATAC peaks.
#'
#' @param RNAseq A numeric matrix containing the RNAseq data
#' @param RNAGroupBy A vector of grouping information for the RNAseq data
#' @param smallRNAseq A numeric matrix containing the smallRNAseq data
#' @param smallRNAGroupBy A vector of grouping information for the smallRNAseq
#'                        data
#' @param proteomics A numeric matrix containing the proteomics data
#' @param proteomicsGroupBy A vector of grouping information for the proteomics
#'                          data
#' @param peakCond1 A data frame containing the differentially accessible ATAC
#'                  peaks for condition 1
#' @param peakCond2 A data frame containing the differentially accessible ATAC
#'                  peaks for condition 2
#'
#' @return NULL
#'
#' @references
#' \insertRef{genomictools}{IntegraTRN}
#'
validateMOInputs <- function(RNAseq,
                             RNAGroupBy,
                             smallRNAseq,
                             smallRNAGroupBy,
                             proteomics,
                             proteomicsGroupBy,
                             peakCond1,
                             peakCond2) {
  # Validate the RNAseq data
  if (!is.null(RNAseq)) {
    validateMatrix(RNAseq, RNAGroupBy)
  } else {
    # Do nothing
  }
  # Validate the smallRNAseq data
  if (!is.null(smallRNAseq)) {
    validateMatrix(smallRNAseq, smallRNAGroupBy)
  } else {
    # Do nothing
  }
  # Validate the proteomics data
  if (!is.null(proteomics)) {
    validateMatrix(proteomics, proteomicsGroupBy)
  } else {
    # Do nothing
  }
  # Replicates are required for the sequencing data
  if ((is.matrix(RNAseq) && ncol(RNAseq) < 2) ||
    (is.matrix(smallRNAseq) && ncol(smallRNAseq) < 2) ||
    (is.matrix(proteomics) && ncol(proteomics) < 2)) {
    stop("At least two replicates are required for each sequencing data.")
  } else {
    # Do nothing
  }
  # Check for missing values in the chromosome information in the ATAC peaks
  if (!is.null(peakCond1) && !is.null(peakCond2)) {
    validatePeaks(peakCond1)
    validatePeaks(peakCond2)
  } else {
    # Do nothing
  }
  return(invisible(NULL))
}

#' Setting non-mRNA omics data to MOList object
#'
#' @keywords internal
#'
#' @description This function sets the non-mRNA omics data to the MOList object.
#'              This is a helper function used by the contructor of the MOList
#'              object to overwrite existing omics data. This is one mode of
#'              the constructor function, which is used to append/exchange
#'              omics data to the MOList object.
#'
#' @param objMOList An object of class MOList for appending/exchanging omics
#'                  data
#' @param smallRNAseq A numeric matrix containing the smallRNAseq data
#' @param smallRNAGroupBy A vector of grouping information for the smallRNAseq
#'                        data
#' @param proteomics A numeric matrix containing the proteomics data
#' @param proteomicsGroupBy A vector of grouping information for the proteomics
#'                          data
#' @param peakCond1 A data frame containing the differentially accessible ATAC
#'                  peaks for condition 1
#' @param peakCond2 A data frame containing the differentially accessible ATAC
#'                  peaks for condition 2
#'
#' @return An object of class MOList
#' \itemize{
#'  \item \code{RNAseq}: A numeric matrix containing the RNAseq data
#'  \item \code{RNAseqSamples}: A list containing the sample names and grouping
#'                            information for the RNAseq data
#'  \item \code{smallRNAseq}: A numeric matrix containing the smallRNAseq data
#'  \item \code{smallRNAseqSamples}: A list containing the sample names and
#'                               grouping information for the smallRNAseq data
#'  \item \code{proteomics}: A numeric matrix containing the proteomics data
#'  \item \code{proteomicsSamples}: A list containing the sample names and
#'                             grouping information for the proteomics data
#'  \item \code{ATACpeaks}: A list containing the differentially accessible ATAC
#'                         peaks for condition 1 and condition 2
#' }
#'
setOmics <- function(objMOList,
                     smallRNAseq,
                     smallRNAGroupBy,
                     proteomics,
                     proteomicsGroupBy,
                     peakCond1,
                     peakCond2) {
  if (!is.null(smallRNAseq)) {
    # Alter the user of overwriting existing data
    if (dim(objMOList@smallRNAseq)[1] != 0) {
      cat("Small RNAseq: Previous small RNAseq data has been overwritten.\n")
      objMOList$DEsmallRNAseq <- NULL # Remove the previous DE results
    } else {
      cat("Small RNAseq: Small RNAseq data has been added.\n")
    }
    objMOList@smallRNAseq <- smallRNAseq
    objMOList@smallRNAseqSamples$samples <- getSampleNames(smallRNAseq)
    objMOList@smallRNAseqSamples$groupBy <- smallRNAGroupBy
  } else {
    cat("Small RNAseq: Optional small RNAseq data has NOT been provided.\n")
  }
  if (!is.null(proteomics)) {
    # Alter again for proteomics data
    if (dim(objMOList@proteomics)[1] != 0) {
      cat("Proteomics:   Previous proteomics data has been overwritten.\n")
      objMOList$DEproteomics <- NULL # Remove the previous DE results
    } else {
      cat("Proteomics:   Proteomics data has been added.\n")
    }
    objMOList@proteomics <- proteomics
    objMOList@proteomicsSamples$samples <- getSampleNames(proteomics)
    objMOList@proteomicsSamples$groupBy <- proteomicsGroupBy
  } else {
    cat("Proteomics:   Optional proteomics data has NOT been provided.\n")
  }
  if (!is.null(peakCond1) && !is.null(peakCond2)) {
    # Same for ATAC peaks
    if (!is.null(objMOList@ATACpeaks$peaksCond1) &&
      !is.null(objMOList@ATACpeaks$peaksCond2)) {
      cat("ATACseq:      Previous ATAC peaks data has been overwritten.\n")
      objMOList$DEATAC <- NULL # Remove the previous DE results
    } else {
      cat("ATACseq:      ATAC peaks data has been added.\n")
    }
    objMOList@ATACpeaks$peaksCond1 <- peakCond1
    objMOList@ATACpeaks$peaksCond2 <- peakCond2
  } else {
    cat("ATACseq:      Optional ATAC peaks data has NOT been provided.\n")
  }
  return(objMOList)
}


#' Construct a new MOList object
#'
#' @keywords internal
#'
#' @description This function constructs a new MOList object, which is used to
#'              store the multi-omics data. Fulfills the CREATION mode of the
#'              constructor function.
#'
#' @param RNAseq A numeric matrix containing the RNAseq data
#' @param RNAGroupBy A vector of grouping information for the RNAseq data. Must
#'                   be a vector of the same length as the number of samples in
#'                   the RNAseq data
#' @param smallRNAseq A numeric matrix containing the smallRNAseq data
#' @param smallRNAGroupBy A vector of grouping information for the smallRNAseq
#'                        data. Must be a vector of the same length as the
#'                        number of samples in the smallRNAseq data
#' @param proteomics A numeric matrix containing the proteomics data
#' @param proteomicsGroupBy A vector of grouping information for the proteomics
#'                          data. Must be a vector of the same length as the
#'                          number of samples in the proteomics data
#' @param peakCond1 A data frame containing the differentially accessible ATAC
#'                  peaks for condition 1. Should follow the BED format
#' @param peakCond2 A data frame containing the differentially accessible ATAC
#'                  peaks for condition 2. Should follow the BED format
#'
#' @return An object of class MOList
#' \itemize{
#' \item \code{RNAseq}: A numeric matrix containing the RNAseq data
#' \item \code{RNAseqSamples}: A list containing the sample names and grouping
#'                           information for the RNAseq data
#' \item \code{smallRNAseq}: A numeric matrix containing the smallRNAseq data
#' \item \code{smallRNAseqSamples}: A list containing the sample names and
#'                             grouping information for the smallRNAseq data
#' \item \code{proteomics}: A numeric matrix containing the proteomics data
#' \item \code{proteomicsSamples}: A list containing the sample names and
#'                            grouping information for the proteomics data
#' \item \code{ATACpeaks}: A list containing the differentially accessible ATAC
#'                         peaks for condition 1 and condition 2
#' }
#'
newMOList <- function(RNAseq,
                      RNAGroupBy,
                      smallRNAseq,
                      smallRNAGroupBy,
                      proteomics,
                      proteomicsGroupBy,
                      peakCond1,
                      peakCond2) {
  # Construct the MOList object
  objMOList <- new("MOList",
    RNAseq = RNAseq,
    RNAseqSamples = list(
      samples = getSampleNames(RNAseq),
      groupBy = RNAGroupBy
    )
  )
  cat("RNAseq:       RNAseq data has been added.\n")
  # Set the non-mRNA omics data
  objMOList <- setOmics(
    objMOList, smallRNAseq, smallRNAGroupBy,
    proteomics, proteomicsGroupBy, peakCond1, peakCond2
  )
  return(objMOList)
}


#' Append/exchange omics data for the existing MOList object
#'
#' @keywords internal
#'
#' @description This function appends/exchanges omics data for the existing
#'              MOList object. Fulfills the APPEND/EXCHANGE mode of the
#'              constructor function.
#'
#' @param objMOList An object of class MOList for appending/exchanging omics
#'                  data
#' @param RNAseq A numeric matrix containing the RNAseq data
#' @param RNAGroupBy A vector of grouping information for the RNAseq data
#' @param smallRNAseq A numeric matrix containing the smallRNAseq data
#' @param smallRNAGroupBy A vector of grouping information for the smallRNAseq
#'                        data
#' @param proteomics A numeric matrix containing the proteomics data
#' @param proteomicsGroupBy A vector of grouping information for the proteomics
#'                          data
#' @param peakCond1 A data frame containing the differentially accessible ATAC
#'                  peaks for condition 1
#' @param peakCond2 A data frame containing the differentially accessible ATAC
#'                  peaks for condition 2
#'
#' @return An object of class MOList
#' \itemize{
#' \item \code{RNAseq}: A numeric matrix containing the RNAseq data
#' \item \code{RNAseqSamples}: A list containing the sample names and grouping
#'                          information for the RNAseq data
#' \item \code{smallRNAseq}: A numeric matrix containing the smallRNAseq data
#' \item \code{smallRNAseqSamples}: A list containing the sample names and
#'                            grouping information for the smallRNAseq data
#' \item \code{proteomics}: A numeric matrix containing the proteomics data
#' \item \code{proteomicsSamples}: A list containing the sample names and
#'                           grouping information for the proteomics data
#' \item \code{ATACpeaks}: A list containing the differentially accessible ATAC
#'                         peaks for condition 1 and condition 2
#' }
#'
modifyMOList <- function(objMOList,
                         RNAseq,
                         RNAGroupBy,
                         smallRNAseq,
                         smallRNAGroupBy,
                         proteomics,
                         proteomicsGroupBy,
                         peakCond1,
                         peakCond2) {
  if (!is.null(RNAseq)) {
    objMOList@RNAseq <- RNAseq
    objMOList@RNAseqSamples$samples <- getSampleNames(RNAseq)
    objMOList@RNAseqSamples$groupBy <- RNAGroupBy
    cat("RNAseq:       Previous RNAseq data has been overwritten.\n")
    objMOList$DERNAseq <- NULL # Remove the previous DE results
  } else {
    # Do nothing
  }
  # Set the non-mRNA omics data
  objMOList <- setOmics(
    objMOList, smallRNAseq, smallRNAGroupBy,
    proteomics, proteomicsGroupBy, peakCond1, peakCond2
  )
  return(objMOList)
}


#' Validate omics raw data with sample information
#'
#' @keywords internal
#'
#' @description This function validates the omics raw data with sample. This is
#'              a helper function used as part of the validator of the MOList
#'              object. Hence, it provides basic sanity checks for the omics
#'              data and sample information that are less stringent than the
#'              input validation checks.
#'
#' @param dataMatrix A numeric matrix containing the omics data
#' @param sampleInfo A list containing the sample names and grouping information
#'                   for the omics data
#'
#' @return NULL
#'
validateOmics <- function(dataMatrix, sampleInfo) {
  if (xor(is.null(dataMatrix), is.null(sampleInfo))) {
    stop("Both the omics data and sample information must be provided.")
  } else if (is.null(dataMatrix) && is.null(sampleInfo)) {
    # Allowed, do nothing
  } else if (!is.numeric(dataMatrix)) {
    stop("The omics data must be a numeric matrix.")
  } else if (ncol(dataMatrix) != length(sampleInfo$samples)) {
    stop("The number of samples in the omics data must match the length of the
    sample information.")
  } else if (any(is.na(dataMatrix))) {
    stop("The omics data contains NA values. Please check the validity of your
    data.")
  } else if (any(is.na(sampleInfo$groupBy)) ||
    ncol(dataMatrix) != length(sampleInfo$groupBy)) {
    stop("Please provide correct grouping information for the omics data.")
  } else {
    # Pass the test, do nothing
  }
  return(invisible(NULL))
}


#' Validator of the MOList S4 object
#'
#' @keywords internal
#'
#' @description This function validates the MOList object. It provides sanity
#'              checks for the omics data and sample information specified in
#'              each slot of the MOList object. All S4 methods of this class
#'              should call this function to validate the MOList object after
#'              any modification is performed.
#'
#' @param objMOList An object of class MOList with any form of modification
#'
#' @return NULL
#'
validateMOList <- function(objMOList) {
  if (class(objMOList)[1] != "MOList") {
    stop("Error validating the MOList class.")
  } else if (is.null(objMOList@RNAseq)) {
    stop("The RNAseq data must be provided.")
  } else {
    # Do nothing
  }
  validateOmics(objMOList@RNAseq, objMOList@RNAseqSamples)
  validateOmics(objMOList@smallRNAseq, objMOList@smallRNAseqSamples)
  validateOmics(objMOList@proteomics, objMOList@proteomicsSamples)
  if (xor(
    is.null(objMOList@ATACpeaks$peaksCond1),
    is.null(objMOList@ATACpeaks$peaksCond2)
  )) {
    stop("Both ATAC peaks files must be provided.")
  } else if (is.null(objMOList@ATACpeaks$peaksCond1) &&
    is.null(objMOList@ATACpeaks$peaksCond2)) {
    # Allowed, do nothing
  } else if (any(is.na(objMOList@ATACpeaks$peaksCond1[, CHROMINFO])) ||
    any(is.na(objMOList@ATACpeaks$peaksCond2[, CHROMINFO]))) {
    stop("Missing values in the chromosome information in the ATAC peaks.")
  } else {
    # Pass the test, do nothing
  }
  return(invisible(NULL))
}


#' Constructor for the MOList object
#'
#' @aliases MOList
#'
#' @description This function is a constructor of the MOList object, which is
#'              used to store the multi-omics data. The MOList object can be
#'              constructed by providing the RNAseq, smallRNAseq,
#'              proteomics, and ATAC peaks data at once, or appending/exchanging
#'              each omics data separately. The RNAseq data is required for the
#'              initial construction of the MOList object. If data of existing
#'              omics type is provided, old data will be replaced.
#'
#' @param objMOList An object of class MOList for appending/exchanging omics
#'                  data. If provided, the MOList object will be modified by
#'                  appending/exchanging the omics data. If not provided, a new
#'                  MOList object will be constructed. Provision of the
#'                  objMOList parameter determines the mode of the constructor
#'                  function
#' @param RNAseq A numeric matrix containing the RNAseq data, with row names as
#'               gene names and column names as sample names
#' @param RNAGroupBy A vector of grouping information for the RNAseq data, used
#'                   to perform differential expression analysis. Must be a
#'                   vector of the same length as the number of samples in the
#'                   RNAseq data
#' @param smallRNAseq A numeric matrix containing the smallRNAseq data, with row
#'                    names as gene names and column names as sample names
#' @param smallRNAGroupBy A vector of grouping information for the smallRNAseq
#'                       data, used to perform differential expression analysis.
#'                        Must be a vector of the same length as the number of
#'                        samples in the smallRNAseq data
#' @param proteomics A numeric matrix containing the proteomics data, with row
#'                   names as gene names and column names as sample names
#' @param proteomicsGroupBy A vector of grouping information for the proteomics
#'                          data, used to perform differential expression
#'                          analysis. Must be a vector of the same length as
#'                          the number of samples in the proteomics data
#' @param ATACpeak1 A character string containing the path to the BED file
#'                  containing the unified ATAC peaks for condition 1. The
#'                  BED file should NOT contain a header line. Alternatively,
#'                  provide a data frame containing the ATAC peaks for condition
#'                  with format consistent with the BED definition
#' @param ATACpeak2 A character string containing the path to the BED file
#'                  containing the unified ATAC peaks for condition 2. The
#'                  BED file should NOT contain a header line. Alternatively,
#'                  provide a data frame containing the ATAC peaks for condition
#'                  with format consistent with the BED definition
#' @note The users should ensure that each BED file used as input
#'                 contains the chromosome regions that are found to have
#'                 increased accessibility in each condition. The BED regions
#'                 can be constructed either by merging the peaks from all
#'                 replicates in each condition using a variable-length method
#'                 or a fixed-width approach, as indicated in Grandi et al. 2022
#'
#' @return An object of class MOList
#' @export MOList
#'
#' @references
#' \insertRef{genomictools}{IntegraTRN}
#'
#' @examples
#' # Generating some example data
#' # Note that the exact value could differ based on the random seed
#' RNAseq <- matrix(sample(1:100, 100, replace = TRUE), ncol = 10)
#' colnames(RNAseq) <- paste0("sample_", seq_len(ncol(RNAseq)))
#' rownames(RNAseq) <- paste0("gene_", seq_len(nrow(RNAseq)))
#' RNAGroupBy <- rep(c("A", "B"), each = 5)
#'
#' smallRNAseq <- matrix(sample(1:100, 20, replace = TRUE), ncol = 4)
#' smallRNAGroupBy <- rep(c("A", "B"), each = 2)
#'
#' proteomics <- matrix(sample(1:100, 30, replace = TRUE), ncol = 6)
#' proteomicsGroupBy <- rep(c("A", "B"), each = 3)
#'
#' # Example 1: Constructing the MOList object at once
#'
#' # Constructing the MOList object at once
#' objMOList1 <- MOList(
#'   RNAseq = RNAseq,
#'   RNAGroupBy = RNAGroupBy,
#'   smallRNAseq = smallRNAseq,
#'   smallRNAGroupBy = smallRNAGroupBy,
#'   proteomics = proteomics,
#'   proteomicsGroupBy = proteomicsGroupBy
#' )
#'
#' # Example 2: Constructing the MOList object step by step
#'
#' # The MOList object can be constructed minimally with the RNAseq data
#' # Further omics data can be appended/exchanged later
#' objMOList2 <- MOList(RNAseq = RNAseq, RNAGroupBy = RNAGroupBy)
#' RNAseq2 <- matrix(sample(1:100, 100, replace = TRUE), ncol = 10)
#' colnames(RNAseq2) <- paste0("sample_", seq_len(ncol(RNAseq2)))
#' objMOList2 <- MOList(objMOList2, RNAseq = RNAseq2, RNAGroupBy = RNAGroupBy)
#' objMOList2 <- MOList(objMOList2,
#'   proteomics = proteomics,
#'   proteomicsGroupBy = proteomicsGroupBy
#' )
#'
#' # Example 3: Exchanging some inputs for the MOList object
#'
#' # The MOList object can be modified by exchanging some inputs
#' objMOList3 <- objMOList2
#' RNAseq3 <- RNAseq
#' RNAseq3[1:10, 1:5] <- 0
#' objMOList3 <- MOList(objMOList3, RNAseq = RNAseq3, RNAGroupBy = RNAGroupBy)
#'
MOList <- function(objMOList = NULL,
                   RNAseq = NULL,
                   RNAGroupBy = NULL,
                   smallRNAseq = NULL,
                   smallRNAGroupBy = NULL,
                   proteomics = NULL,
                   proteomicsGroupBy = NULL,
                   ATACpeak1 = NULL,
                   ATACpeak2 = NULL) {
  # One of objMOList and RNAseq must be given, which dictates whether the
  # constructor creates a new object or appends/exchanges omics data
  if (is.null(objMOList) && is.null(RNAseq)) {
    stop("Please provide either an object of class MOList or the RNAseq data.")
  } else {
    # Do nothing
  }

  # Read the ATAC peaks data if provided and are valid paths
  if (is.data.frame(ATACpeak1) && is.data.frame(ATACpeak2)) {
    peakCond1 <- ATACpeak1
    colnames(peakCond1)[1:3] <- CHROMINFO
    peakCond2 <- ATACpeak2
    colnames(peakCond2)[1:3] <- CHROMINFO
  } else if (!is.null(ATACpeak1) && !is.null(ATACpeak2) &&
    is.character(ATACpeak1) && is.character(ATACpeak2)) {
    if (file.exists(ATACpeak1) && file.exists(ATACpeak2)) {
      peakCond1 <- utils::read.table(ATACpeak1, header = FALSE, sep = "\t")
      colnames(peakCond1)[1:3] <- CHROMINFO
      peakCond2 <- utils::read.table(ATACpeak2, header = FALSE, sep = "\t")
      colnames(peakCond2)[1:3] <- CHROMINFO
    } else {
      stop("Please provide valid paths to the ATAC peaks files.")
    }
  } else if (xor(is.null(ATACpeak1), is.null(ATACpeak2))) {
    stop("Please provide both valid ATAC peaks files.")
  } else {
    peakCond1 <- NULL
    peakCond2 <- NULL
  }

  # Input validation
  validateMOInputs(
    RNAseq, RNAGroupBy, smallRNAseq, smallRNAGroupBy,
    proteomics, proteomicsGroupBy, peakCond1, peakCond2
  )

  # Construct or modify the MOList object
  if (is.null(objMOList)) {
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
  # Final validation of the MOList object
  validateMOList(newObjMOList)
  return(newObjMOList)
}


# Define a set of setters for the data slots of the MOList object
# Must be used in caution, the user should not have access to these functions

#' Setter for the RNAseq slot
#'
#' @aliases RNAseq<-
#'
#' @keywords internal
#'
#' @description This function is a setter for the RNAseq slot of the MOList.
#'              The user are NOT allowed to use this function directly, since
#'              use without caution may cause the MOList object to be invalid.
#'
#' @param x An object of class MOList for appending/exchanging omics data
#' @param value A numeric matrix containing the RNAseq data
#'
#' @return An object of class MOList
#' \itemize{
#' \item \code{RNAseq}: A numeric matrix containing the RNAseq data
#' \item \code{RNAseqSamples}: A list containing the sample names and grouping
#'                         information for the RNAseq data
#' \item \code{smallRNAseq}: A numeric matrix containing the smallRNAseq data
#' \item \code{smallRNAseqSamples}: A list containing the sample names and
#'                              grouping information for the smallRNAseq data
#' \item \code{proteomics}: A numeric matrix containing the proteomics data
#' \item \code{proteomicsSamples}: A list containing the sample names and
#'                           grouping information for the proteomics data
#' \item \code{ATACpeaks}: A list containing the differentially accessible ATAC
#'                         peaks for condition 1 and condition 2
#' }
#'
setGeneric("RNAseq<-", function(x, value) standardGeneric("RNAseq<-"))
setMethod("RNAseq<-", "MOList", function(x, value) {
  if (is.null(value) || any(colnames(value) != x@RNAseqSamples$samples)) {
    stop("The sample names of the smallRNAseq data must match the RNAseq data.")
  } else {
    x@RNAseq <- value
    validateMOList(x)
    return(x)
  }
})

#' Setter for the smallRNAseq slot
#'
#' @aliases smallRNAseq<-
#'
#' @keywords internal
#'
#' @param x An object of class MOList for appending/exchanging omics data
#' @param value A numeric matrix containing the smallRNAseq data
#'
#' @return An object of class MOList
#' \itemize{
#' \item \code{RNAseq}: A numeric matrix containing the RNAseq data
#' \item \code{RNAseqSamples}: A list containing the sample names and grouping
#'                         information for the RNAseq data
#' \item \code{smallRNAseq}: A numeric matrix containing the smallRNAseq data
#' \item \code{smallRNAseqSamples}: A list containing the sample names and
#'                              grouping information for the smallRNAseq data
#' \item \code{proteomics}: A numeric matrix containing the proteomics data
#' \item \code{proteomicsSamples}: A list containing the sample names and
#'                           grouping information for the proteomics data
#' \item \code{ATACpeaks}: A list containing the differentially accessible ATAC
#'                         peaks for condition 1 and condition 2
#' }
#'
setGeneric(
  "smallRNAseq<-",
  function(x, value) standardGeneric("smallRNAseq<-")
)
setMethod("smallRNAseq<-", "MOList", function(x, value) {
  if (is.null(value) || any(colnames(value) != x@smallRNAseqSamples$samples)) {
    stop("The sample names of the smallRNAseq data must match the RNAseq data.")
  } else {
    x@smallRNAseq <- value
    validateMOList(x)
    return(x)
  }
})

#' Setter for the proteomics slot
#'
#' @aliases proteomics<-
#'
#' @keywords internal
#'
#' @param x An object of class MOList for appending/exchanging omics data
#' @param value A numeric matrix containing the proteomics data
#'
#' @return An object of class MOList
#' \itemize{
#' \item \code{RNAseq}: A numeric matrix containing the RNAseq data
#' \item \code{RNAseqSamples}: A list containing the sample names and grouping
#'                         information for the RNAseq data
#' \item \code{smallRNAseq}: A numeric matrix containing the smallRNAseq data
#' \item \code{smallRNAseqSamples}: A list containing the sample names and
#'                              grouping information for the smallRNAseq data
#' \item \code{proteomics}: A numeric matrix containing the proteomics data
#' \item \code{proteomicsSamples}: A list containing the sample names and
#'                           grouping information for the proteomics data
#' \item \code{ATACpeaks}: A list containing the differentially accessible ATAC
#'                         peaks for condition 1 and condition 2
#' }
#'
setGeneric(
  "proteomics<-",
  function(x, value) standardGeneric("proteomics<-")
)
setMethod("proteomics<-", "MOList", function(x, value) {
  if (is.null(value) || any(colnames(value) != x@proteomicsSamples$samples)) {
    stop("The sample names of the proteomics data must match the RNAseq data.")
  } else {
    x@proteomics <- value
    validateMOList(x)
    return(x)
  }
})

#' Setter for the ATACpeaks slot
#'
#' @aliases ATACpeaks<-
#'
#' @keywords internal
#'
#' @param x An object of class MOList for appending/exchanging omics data
#' @param value A list containing the differentially accessible ATAC peaks for
#'              condition 1 and condition 2
#'
#' @return An object of class MOList
#' \itemize{
#' \item \code{RNAseq}: A numeric matrix containing the RNAseq data
#' \item \code{RNAseqSamples}: A list containing the sample names and grouping
#'                         information for the RNAseq data
#' \item \code{smallRNAseq}: A numeric matrix containing the smallRNAseq data
#' \item \code{smallRNAseqSamples}: A list containing the sample names and
#'                              grouping information for the smallRNAseq data
#' \item \code{proteomics}: A numeric matrix containing the proteomics data
#' \item \code{proteomicsSamples}: A list containing the sample names and
#'                           grouping information for the proteomics data
#' \item \code{ATACpeaks}: A list containing the differentially ATAC peaks for
#'                         condition 1 and condition 2
#' }
#'
setGeneric(
  "ATACpeaks<-",
  function(x, value) standardGeneric("ATACpeaks<-")
)
setMethod("ATACpeaks<-", "MOList", function(x, value) {
  if (is.null(value)) {
    stop("Invalid ATAC peak list.")
  } else {
    x@ATACpeaks <- value
    validateMOList(x)
    return(x)
  }
})


# Define a set of getters for the data slots of the MOList object, allowing
# the user to retrieve information from the MOList object

#' Getter for the count-based omics data from the MOList object
#'
#' @aliases getRawData
#'
#' @aliases getRawData,MOList-method
#'
#' @description This function is a getter for the count-based omics data from
#'              the MOList object. The user can retrieve the RNAseq,
#'              smallRNAseq, and proteomics data from the MOList object.
#'
#' @param x An object of class MOList for retrieving omics data
#' @param omics A character string specifying the omics data to be retrieved
#'        must be one of "RNAseq", "smallRNAseq", "proteomics", or "ATAC"
#'
#' @return Raw data of the specified omics type: a numeric matrix for
#'         count-based omics data, and a list of data frames for ATAC peaks
#'         that are differentially accessible between two conditions
#'
#' @export
#'
#' @examples
#' # Use the package-provided example MOList object
#' data("expMOList")
#'
#' # Example 1: Retrieving the RNAseq raw counts
#' getRawData(expMOList, omics = "RNAseq")
#'
#' # Example 2: Retrieving the smallRNAseq raw counts
#' getRawData(expMOList, omics = "smallRNAseq")
#'
#' # Example 3: Retrieving the proteomics raw counts
#' getRawData(expMOList, omics = "proteomics")
#'
#' # Example 4: Retrieving the ATAC peaks as a list
#' getRawData(expMOList, omics = "ATAC")
#'
setGeneric(
  "getRawData",
  function(x, omics) standardGeneric("getRawData")
)
setMethod("getRawData", "MOList", function(x, omics) {
  omicData <- switch(omics,
    RNAseq = x@RNAseq,
    smallRNAseq = x@smallRNAseq,
    proteomics = x@proteomics,
    ATAC = x@ATACpeaks
  )
  return(omicData)
})


#' Retrieving the sample information from the MOList object
#'
#' @keywords internal
#'
#' @param x A MOList object containing the omics data
#' @param experiment A character string specifying the experiment type, must
#'                   be one of "RNAseq", "smallRNAseq", and "proteomics"
#'
#' @return A list containing the sample names and grouping information for the
#'         specified omics data
#'
setGeneric(
  "getSampleInfo",
  function(x, experiment) standardGeneric("getSampleInfo")
)
setMethod("getSampleInfo", "MOList", function(x, experiment) {
  sampleInfo <- switch(experiment,
    RNAseq = x@RNAseqSamples,
    smallRNAseq = x@smallRNAseqSamples,
    proteomics = x@proteomicsSamples
  )
  return(sampleInfo)
})


#' print MOList object
#'
#' @aliases print,MOList-method
#'
#' @description This function prints the MOList object
#'
#' @param object An object of the DETag class
#'
#' @return NULL
#'
#' @export
#'
#' @examples
#' # Use the package-provided example MOList object
#' data("expMOList")
#'
#' # Print the MOList object
#' print(expMOList)
#'
#' # or simply type the object name
#' expMOList
#'
setMethod("show", "MOList", function(object) {
  cat("MOList object with the following slots:\n")
  cat(
    "RNAseq with", nrow(object@RNAseq), "genes and",
    ncol(object@RNAseq), "samples\n"
  )
  print(utils::head(object@RNAseq))
  cat("\n")

  cat(
    "Small RNAseq with", nrow(object@smallRNAseq), "genes and",
    ncol(object@smallRNAseq), "samples\n"
  )
  print(utils::head(object@smallRNAseq))
  cat("\n")

  cat(
    "Proteomics with", nrow(object@proteomics), "genes and",
    ncol(object@proteomics), "samples\n"
  )
  print(utils::head(object@proteomics))
  cat("\n")

  cat("ATACseq:\n")
  cat("Peaks for condition 1:", nrow(object@ATACpeaks$peaksCond1), "\n")
  print(utils::head(object@ATACpeaks$peaksCond1))
  cat("\n")
  cat("Peaks for condition 2:", nrow(object@ATACpeaks$peaksCond2), "\n")
  print(utils::head(object@ATACpeaks$peaksCond2))
  cat("\n")

  return(invisible(NULL))
})


#' Setting conversion between protein and gene names
#'
#' @aliases setGene2Protein,MOList-method
#'
#' @description This function sets the conversion between protein and gene names
#'              for the proteomics data. The user can provide a data frame
#'              containing the conversion information. The data frame must
#'              contain two columns, one for the protein names and the other
#'              for the gene names. This must be provided to effectively utilize
#'              the proteomics data. The users are responsible for the validity
#'              of the conversion information.
#'
#' @param x A MOList object containing the omics data
#' @param conversion A data frame containing the conversion information between
#'                   protein and gene names, with names matching to the
#'                   proteomics  and RNAseq raw count data
#' \itemize{
#' \item \code{protein}: Name of the protein
#' \item \code{gene}: Name of the gene
#' }
#'
#' @return An object of class MOList
#'
#' @export
#'
#' @examples
#' # Create example RNAseq and proteomics data
#' rnaseq <- matrix(sample(1:100, 50, replace = TRUE), ncol = 10)
#' rownames(rnaseq) <- paste0("gene_", seq_len(nrow(rnaseq)))
#' rnaGroupBy <- rep(c("A", "B"), each = 5)
#' proteomics <- matrix(sample(1:100, 30, replace = TRUE), ncol = 6)
#' rownames(proteomics) <- paste0("protein_", seq_len(nrow(proteomics)))
#' proteomicsGroupBy <- rep(c("A", "B"), each = 3)
#'
#' # Create an example MOList object
#' objMOList <- MOList(
#'   RNAseq = rnaseq,
#'   RNAGroupBy = rnaGroupBy,
#'   proteomics = proteomics,
#'   proteomicsGroupBy = proteomicsGroupBy
#' )
#'
#' # Create an example conversion information
#' conversion <- data.frame(
#'   protein = paste0("protein_", seq_len(nrow(rnaseq))),
#'   gene = paste0("gene_", seq_len(nrow(proteomics)))
#' )
#'
#' # Set the conversion information to the MOList object
#' objMOList <- setGene2Protein(objMOList, conversion)
#'
setGeneric(
  "setGene2Protein",
  function(x, conversion) standardGeneric("setGene2Protein")
)
setMethod("setGene2Protein", "MOList", function(x, conversion) {
  if (is.null(x@proteomics)) {
    stop("Please provide the proteomics data.")
  } else if (!is.data.frame(conversion)) {
    stop("The conversion information must be a data frame.")
  } else if (!all(c("protein", "gene") %in% colnames(conversion))) {
    stop("The conversion information must contain the protein and gene names.")
  } else {
    # Continue
  }
  x$gene2protein <- conversion
  return(x)
})


# [END]
