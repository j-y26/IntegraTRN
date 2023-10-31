# Purpose: Creating the Multi-Omics List S4 class
# Author: Jielin Yang
# Date: 2023-10-29
# Version: 1.0
# Bugs and Issues: None


# Global variables
CHROMINFO <- c("chrom", "chromStart", "chromEnd")


#' Multi-Omics List (MOList) S4 class
#'
#' @description This class is used to store the multi-omics data, including
#'              RNAseq, smallRNAseq, proteomics, and ATAC peaks data.
#' @slot RNAseq A numeric matrix containing the RNAseq data
#' @slot RNAseqSamples A list containing the sample names and grouping
#'                     information for the RNAseq data
#' @slot smallRNAseq A numeric matrix containing the smallRNAseq data
#' @slot smallRNAseqSamples A list containing the sample names and grouping
#'                          information for the smallRNAseq data
#' @slot proteomics A numeric matrix containing the proteomics data
#' @slot proteomicsSamples A list containing the sample names and grouping
#'                         information for the proteomics data
#' @slot ATACpeaks A list containing the ATAC peaks for condition 1 and
#'                 condition 2
#' @slot .Data A list containing differential expression analysis data, used
#'             in the same way as a list
#'
#' @importFrom methods setClass
#'
#' @exportClass MOList
#'
#' @references
#' Advanced R by H. Wickham. Access: https://adv-r.hadley.nz/index.html
#'
methods::setClass("MOList",
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
#'              the number of samples in the sequencing data.
#' @param matrix A numeric matrix containing the sequencing data
#' @param groupBy A vector of grouping information for the sequencing data
#'
#' @return NULL
#'
#' @references
#' Yang Liao, Gordon K. Smyth, Wei Shi, featureCounts: an efficient general
#' purpose program for assigning sequence reads to genomic features,
#' Bioinformatics, Volume 30, Issue 7, April 2014, Pages 923–930,
#' https://doi.org/10.1093/bioinformatics/btt656
#'
#' @examples
#' \dontrun{
#' # Assuming RNAseq, smallRNAseq, and proteomics are matrices of count data
#'
#' # Validating the RNAseq data
#' validateMatrix(RNAseq, RNAGroupBy)
#' # Validating the smallRNAseq data
#' validateMatrix(smallRNAseq, smallRNAGroupBy)
#' # Validating the proteomics data
#' validateMatrix(proteomics, proteomicsGroupBy)
#' }
#'
validateMatrix <- function(matrix, groupBy) {
  if (any(is.na(matrix))) {
    stop("The sequencing data contains NA values. Please check the validity of
    your data.")
  } else if (any(is.na(groupBy))) {
    stop("Please provide correct grouping information for the sequencing data.")
  } else if (!is.numeric(matrix) || !is.matrix(matrix)) {
    stop("The sequencing data must be a numeric matrix.")
  } else if (ncol(matrix) != length(groupBy)) {
    stop("The number of samples in the sequencing data must match the length of
    the grouping information.")
  } else {
    # Do nothing
  }
  return(invisible(NULL))
}


#' Validate the inputs for the MOList constructor
#'
#' @keywords internal
#'
#' @description This function validates the inputs for the MOList constructor.
#'              It uses the validateMatrix function to validate the sequencing
#'              data, and checks for missing values in the chromosome
#'              information in the ATAC peaks.
#' @param RNAseq A numeric matrix containing the RNAseq data
#' @param RNAGroupBy A vector of grouping information for the RNAseq data
#' @param smallRNAseq A numeric matrix containing the smallRNAseq data
#' @param smallRNAGroupBy A vector of grouping information for the smallRNAseq
#'                        data
#' @param proteomics A numeric matrix containing the proteomics data
#' @param proteomicsGroupBy A vector of grouping information for the proteomics
#'                          data
#' @param peakCond1 A data frame containing the ATAC peaks for condition 1
#' @param peakCond2 A data frame containing the ATAC peaks for condition 2
#'
#' @return NULL
#'
#' @references
#' Fischer, D. (2020). GenomicTools.fileHandler: File handlers for genomic data
#' analysis.
#'
#' @examples
#' \dontrun{
#' # Assuming RNAseq, smallRNAseq, and proteomics are matrices of count data
#'
#' # Validating the inputs for the MOList constructor
#' validateMOInputs(
#'   RNAseq, RNAGroupBy, smallRNAseq, smallRNAGroupBy,
#'   proteomics, proteomicsGroupBy, peakCond1, peakCond2
#' )
#' }
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
    if (any(is.na(peakCond1[, CHROMINFO])) ||
      any(is.na(peakCond2[, CHROMINFO]))) {
      stop("Missing values in the chromosome information in the ATAC peaks.")
    } else {
      # Do nothing
    }
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
#'
#' @param objMOList An object of class MOList for appending/exchanging omics
#'                  data
#' @param smallRNAseq A numeric matrix containing the smallRNAseq data
#' @param smallRNAGroupBy A vector of grouping information for the smallRNAseq
#'                        data
#' @param proteomics A numeric matrix containing the proteomics data
#' @param proteomicsGroupBy A vector of grouping information for the proteomics
#'                          data
#' @param peakCond1 A data frame containing the ATAC peaks for condition 1
#' @param peakCond2 A data frame containing the ATAC peaks for condition 2
#'
#' @return An object of class MOList
#' \itemize{
#'  \item \code{RNAseq}: A numeric matrix containing the RNAseq data
#'  \item \code{RNAseqSamples}: A list containing the sample names and grouping
#'                            information for the RNAseq data
#'  \item \code{smallRNAseq}: A numeric matrix containing the smallRNAseq data
#' \item \code{smallRNAseqSamples}: A list containing the sample names and
#'                               grouping information for the smallRNAseq data
#' \item \code{proteomics}: A numeric matrix containing the proteomics data
#' \item \code{proteomicsSamples}: A list containing the sample names and
#'                             grouping information for the proteomics data
#' \item \code{ATACpeaks}: A list containing the ATAC peaks for condition 1 and
#'                      condition 2
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming RNAseq, smallRNAseq, and proteomics are matrices of count data
#'
#' # Setting the omics data to their respective slot in the MOList object
#' objMOList <- setOmics(
#'   objMOList, smallRNAseq, smallRNAGroupBy,
#'   proteomics, proteomicsGroupBy, peakCond1, peakCond2
#' )
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
    objMOList@smallRNAseq <- smallRNAseq
    objMOList@smallRNAseqSamples$samples <- getSampleNames(smallRNAseq)
    objMOList@smallRNAseqSamples$groupBy <- smallRNAGroupBy
  } else {
    # Do nothing
  }
  if (!is.null(proteomics)) {
    objMOList@proteomics <- proteomics
    objMOList@proteomicsSamples$samples <- getSampleNames(proteomics)
    objMOList@proteomicsSamples$groupBy <- proteomicsGroupBy
  } else {
    # Do nothing
  }
  if (!is.null(peakCond1) && !is.null(peakCond2)) {
    objMOList@ATACpeaks$peaksCond1 <- peakCond1
    objMOList@ATACpeaks$peaksCond2 <- peakCond2
  } else {
    # Do nothing
  }
  return(objMOList)
}


#' Construct a new MOList object
#'
#' @keywords internal
#'
#' @param RNAseq A numeric matrix containing the RNAseq data
#' @param RNAGroupBy A vector of grouping information for the RNAseq data
#' @param smallRNAseq A numeric matrix containing the smallRNAseq data
#' @param smallRNAGroupBy A vector of grouping information for the smallRNAseq
#'                        data
#' @param proteomics A numeric matrix containing the proteomics data
#' @param proteomicsGroupBy A vector of grouping information for the proteomics
#'                          data
#' @param peakCond1 A data frame containing the ATAC peaks for condition 1
#' @param peakCond2 A data frame containing the ATAC peaks for condition 2
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
#' \item \code{ATACpeaks}: A list containing the ATAC peaks for condition 1 and
#'                     condition 2
#' }
#'
#' @importFrom methods new
#'
#' @references
#' Advanced R by H. Wickham. Access: https://adv-r.hadley.nz/index.html
#'
#' @examples
#' \dontrun{
#' # Assuming RNAseq, smallRNAseq, and proteomics are matrices of count data
#' # and RNAGroupBy, smallRNAGroupBy, and proteomicsGroupBy are vectors of
#' # grouping information for the sequencing data
#'
#' # Constructing a new MOList object
#' objMOList <- newMOList(
#'   RNAseq, RNAGroupBy, smallRNAseq, smallRNAGroupBy,
#'   proteomics, proteomicsGroupBy, peakCond1, peakCond2
#' )
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
  objMOList <- methods::new("MOList",
    RNAseq = RNAseq,
    RNAseqSamples = list(
      samples = getSampleNames(RNAseq),
      groupBy = RNAGroupBy
    )
  )
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
#' @param objMOList An object of class MOList for appending/exchanging omics data
#' @param RNAseq A numeric matrix containing the RNAseq data
#' @param RNAGroupBy A vector of grouping information for the RNAseq data
#' @param smallRNAseq A numeric matrix containing the smallRNAseq data
#' @param smallRNAGroupBy A vector of grouping information for the smallRNAseq
#'                        data
#' @param proteomics A numeric matrix containing the proteomics data
#' @param proteomicsGroupBy A vector of grouping information for the proteomics
#'                          data
#' @param peakCond1 A data frame containing the ATAC peaks for condition 1
#' @param peakCond2 A data frame containing the ATAC peaks for condition 2
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
#' \item \code{ATACpeaks}: A list containing the ATAC peaks for condition 1 and
#'                    condition 2
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming RNAseq, smallRNAseq, and proteomics are matrices of count data
#' # and RNAGroupBy, smallRNAGroupBy, and proteomicsGroupBy are vectors of
#' # grouping information for the sequencing data
#'
#' # Modifying an existing MOList object
#' objMOList <- modifyMOList(
#'   objMOList, RNAseq, RNAGroupBy, smallRNAseq,
#'   smallRNAGroupBy, proteomics, proteomicsGroupBy,
#'   peakCond1, peakCond2
#' )
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
#' @param dataMatrix A numeric matrix containing the omics data
#' @param sampleInfo A list containing the sample names and grouping information
#'                   for the omics data
#'
#' @return NULL
#'
#' @examples
#' \dontrun{
#' # Assuming RNAseq is a matrix of count data and RNASamples is a list of
#' # sample names and grouping information for the RNAseq data
#'
#' # Validating the input data
#' validateOmics(RNAseq, RNAsamples)
#' }
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
#' @param objMOList An object of class MOList for appending/exchanging omics data
#'
#' @return NULL
#'
#' @examples
#' \dontrun{
#' # Assuming objMOList is an object of class MOList
#'
#' # Validating the MOList object
#' validateMOList(objMOList)
#' }
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
#' @importFrom GenomicTools.fileHandler importBed
#'
#' @references
#' Fischer, D. (2020). GenomicTools.fileHandler: File handlers for genomic data
#' analysis.
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
                   pathATACpeak1 = NULL,
                   pathATACpeak2 = NULL) {
  # One of objMOList and RNAseq must be given, which dictates whether the
  # constructor creates a new object or appends/exchanges omics data
  if (is.null(objMOList) && is.null(RNAseq)) {
    stop("Please provide either an object of class MOList or the RNAseq data.")
  } else {
    # Do nothing
  }

  # Read the ATAC peaks data if provided
  if (!is.null(pathATACpeak1) && !is.null(pathATACpeak2)) {
    peakCond1 <- GenomicTools.fileHandler::importBed(pathATACpeak1)
    peakCond2 <- GenomicTools.fileHandler::importBed(pathATACpeak2)
  } else if (xor(is.null(pathATACpeak1), is.null(pathATACpeak2))) {
    stop("Please provide both ATAC peaks files.")
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
#' \item \code{ATACpeaks}: A list containing the ATAC peaks for condition 1 and
#'                   condition 2
#' }
#'
#' @importFrom methods setGeneric setMethod
#'
#' @examples
#' \dontrun{
#' # Assuming RNAseq is a matrix of count data
#'
#' # Setting the RNAseq data to the MOList object
#' RNAseq(myMOList) <- RNAseq
#' }
#'
methods::setGeneric("RNAseq<-", function(x, value) standardGeneric("RNAseq<-"))
methods::setMethod("RNAseq<-", "MOList", function(x, value) {
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
#' \item \code{ATACpeaks}: A list containing the ATAC peaks for condition 1 and
#'                   condition 2
#' }
#'
#' @importFrom methods setGeneric setMethod
#'
#' @examples
#' \dontrun{
#' # Assuming smallRNAseq is a matrix of count data
#'
#' # Setting the smallRNAseq data to the MOList object
#' smallRNAseq(myMOList) <- smallRNAseq
#' }
#'
methods::setGeneric(
  "smallRNAseq<-",
  function(x, value) standardGeneric("smallRNAseq<-")
)
methods::setMethod("smallRNAseq<-", "MOList", function(x, value) {
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
#' \item \code{ATACpeaks}: A list containing the ATAC peaks for condition 1 and
#'                   condition 2
#' }
#'
#' @importFrom methods setGeneric setMethod
#'
#' @examples
#' \dontrun{
#' # Assuming proteomics is a matrix of count data
#'
#' # Setting the proteomics data to the MOList object
#' proteomics(myMOList) <- proteomics
#' }
#'
methods::setGeneric(
  "proteomics<-",
  function(x, value) standardGeneric("proteomics<-")
)
methods::setMethod("proteomics<-", "MOList", function(x, value) {
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
#' @param value A list containing the ATAC peaks for condition 1 and condition 2
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
#' \item \code{ATACpeaks}: A list containing the ATAC peaks for condition 1 and
#'                   condition 2
#' }
#'
#' @importFrom methods setGeneric setMethod
#'
#' @examples
#' \dontrun{
#' # Assuming peakCond1 and peakCond2 are data frames containing the ATAC peaks
#' # for condition 1 and condition 2, respectively
#'
#' # Setting the ATAC peaks data to the MOList object
#' ATACpeaks(myMOList) <- list(peaksCond1 = peakCond1, peaksCond2 = peakCond2)
#' }
#'
methods::setGeneric(
  "ATACpeaks<-",
  function(x, value) standardGeneric("ATACpeaks<-")
)
methods::setMethod("ATACpeaks<-", "MOList", function(x, value) {
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
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Using the example MOList object
#' dataRNAseq <- getCounts(myMOList, "RNAseq")
#' dataSmallRNAseq <- getCounts(myMOList, "smallRNAseq")
#' dataProteomics <- getCounts(myMOList, "proteomics")
#' dataATAC <- getCounts(myMOList, "ATAC")
#' }
#'
methods::setGeneric(
  "getRawData",
  function(x, omics) standardGeneric("getRawData")
)
methods::setMethod("getRawData", "MOList", function(x, omics) {
  omicData <- switch(omics,
    RNAseq = x@RNAseq,
    smallRNAseq = x@smallRNAseq,
    proteomics = x@proteomics,
    ATAC = x@ATACpeaks
  )
  return(omicData)
})


#' Retrieving the grouping information for a specific omics data
#'
#' @keywords internal
#'
#' @param x A MOList object containing the omics data
#' @param experiment A character string specifying the experiment type, must
#'                   be one of "RNAseq", "smallRNAseq", and "proteomics"
#'
#' @return A vector of grouping information
#'
#' @importFrom methods setGeneric setMethod
#'
#' @examples
#' \dontrun{
#' # Assuming myMOList is a MOList object
#'
#' # Get the grouping information for the RNAseq data
#' getSampleInfo(myMOList, "RNAseq")
#'
#' # Get the grouping information for the smallRNAseq data
#' getSampleInfo(myMOList, "smallRNAseq")
#'
#' # Get the grouping information for the proteomics data
#' getSampleInfo(myMOList, "proteomics")
#' }
#'
methods::setGeneric(
  "getSampleInfo",
  function(x, experiment) standardGeneric("getSampleInfo")
)
methods::setMethod("getSampleInfo", "MOList", function(x, experiment) {
  sampleInfo <- switch(experiment,
    RNAseq = x@RNAseqSamples$groupBy,
    smallRNAseq = x@smallRNAseqSamples$groupBy,
    proteomics = x@proteomicsSamples$groupBy
  )
  return(sampleInfo)
})
