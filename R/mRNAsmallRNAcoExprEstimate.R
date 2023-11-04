# Purpose: Predictive estimation of regulatory small RNA - mRNA interaction
#          based on matched co-expression data
# Author: Jielin Yang
# Date: 2023-11-30
# Version: 1.0
# Bugs and Issues: None


#' Validate RNAseq and small RNAseq sample data frames for matching
#' 
#' @description The data frames should contain the same number of samples
#'              and samples names as already documented in the MOList object.
#'              However, they do not need to be in the same order.
#' 
#' @param objMOList An object of class MOList
#' @param sampleDFRNAseq A data frame containing the RNAseq sample information.
#'                       Each sample is a row in the data frame, with row names
#'                       being the sample names that correspond to the sample
#'                       names in the MOList object.
#' @param sampleDFSmallRNAseq A data frame containing the small RNAseq sample
#'                            information. Each sample is a row in the data
#'                            frame, with row names being the sample names that
#'                            correspond to the sample names in the MOList
#'                            object.
#' 
#' @return NULL
#' 
validateSampleDFs <- function(objMOList, sampleDFRNAseq, sampleDFSmallRNAseq) {
  # Check of sample numbers
  if (length(objMOList@RNAseqSamples$samples) != nrow(sampleDFRNAseq) ||
      length(objMOList@SmallRNAseqSamples$samples) != 
                                                    nrow(sampleDFSmallRNAseq)) {
      stop("The number of samples in the sample data frames does not match the 
      number of samples in the MOList object.")
  } else {
    # Continue
  }
  # Check for sample name matches
  if (!all(objMOList@RNAseqSamples$samples %in% rownames(sampleDFRNAseq)) ||
      !all(objMOList@SmallRNAseqSamples$samples %in% 
                                              rownames(sampleDFSmallRNAseq))) {
    stop("The sample names in the sample data frames do not match the sample 
    names in the MOList object.")
  } else {
    # Continue
  }
  return(invisible(NULL))
}


#' Performing nearest neighbor matching for RNAseq and small RNAseq samples
#' 
#' @description This function performs nearest neighbor matching for RNAseq and
#'              small RNAseq samples using two sample data frames. The two data
#'              frames should contain the same number of variables, including
#'              the main grouping variable for DE analysis.
#' 
#' @param sampleDFRNAseq A data frame containing the RNAseq sample information.
#'                       Each sample is a row in the data frame, with row names
#'                       being the sample names that correspond to the sample
#'                       names in the MOList object.
#' @param sampleDFSmallRNAseq A data frame containing the small RNAseq sample
#'                            information. Each sample is a row in the data
#'                            frame, with row names being the sample names that
#'                            correspond to the sample names in the MOList
#'                            object.
#' 
#' @return A list of two numeric vectors with equal length, where each element
#'         in the first vector is the index of the RNAseq sample, and the
#'         corresponding element in the second vector is the index of the
#'         matched small RNAseq sample.
#' 
nnRNAMatch <- function(sampleDFRNAseq, sampleDFSmallRNAseq)



#' Match RNAseq and small RNAseq samples
#' 
#' @description This function generates a one-to-one matching between the
#'              RNAseq samples and the small RNAseq samples using the
#'              nearest neighbor matching algorithm.
#' 
#' @param objMOList An object of class MOList
#' @param sampleDFRNAseq A data frame containing the RNAseq sample information.
#'                       Each sample is a row in the data frame, with row names
#'                       being the sample names that correspond to the sample
#'                       names in the MOList object.
#' @param sampleDFSmallRNAseq A data frame containing the small RNAseq sample
#'                            information. Each sample is a row in the data
#'                            frame, with row names being the sample names that
#'                            correspond to the sample names in the MOList
#'                            object.
#' @param varMatch A vector of variable names that will be used for matching.
#'                 The variable names should be present in both sampleDFRNAseq
#'                 and sampleDFSmallRNAseq. If NULL, default to use all
#'                 common variables in the sample data frames.
#' 
#' @return An object of class MOList with the matching information stored in
#'         a list of two numeric vectors with equal length, where each element
#'         in the first vector is the index of the RNAseq sample, and the
#'         corresponding element in the second vector is the index of the
#'         matched small RNAseq sample. These indices match to the sample names
#'         in the RNAseqSamples and SmallRNAseqSamples slots of the MOList
#'         object.
#' 
matchSamplesRNAsmallRNA <- function(objMOList, 
                         sampleDFRNAseq = NULL, 
                         sampleDFSmallRNAseq = NULL, 
                         varMatch = NULL) {
  # Check if the the MOList object has smallRNAseq data
  if (is.null(getRawData(objMOList, SMALLRNA))) {
    stop("The MOList object does not contain smallRNAseq data.")
  }

  groupByRNA <- objMOList@RNAseqSamples$groupBy
  groupBySmallRNA <- objMOList@SmallRNAseqSamples$groupBy

  # Only use sample information when both sampleDFRNAseq and
  # sampleDFSmallRNAseq are provided
  if (!is.null(sampleDFRNAseq) && !is.null(sampleDFSmallRNAseq)) {
    validateSampleDFs(objMOList, sampleDFRNAseq, sampleDFSmallRNAseq)
    if (!is.null(varMatch)) {
      invalidVars <- setdiff(varMatch, intersect(colnames(sampleDFRNAseq), 
                                                 colnames(sampleDFSmallRNAseq)))
      if (length(invalidVars) > 0) {
        warning(paste("The following variables are not present in both sample 
                      data frames and are ignored:", invalidVars))
        varMatch <- setdiff(varMatch, invalidVars)
      } else {
        # Continue
      }
    } else {
      varMatch <- intersect(colnames(sampleDFRNAseq), 
                            colnames(sampleDFSmallRNAseq))
    }
    # Keep only the usable sample information
    sampleDFRNAseq <- sampleDFRNAseq[, varMatch]
    sampleDFSmallRNAseq <- sampleDFSmallRNAseq[, varMatch]

    # Since the key grouping variable may or may not be present in the selected
    # matching variables, we need to check if it is present, and if not, add it
    # to the matching variables (and the sample data frames)
    # This is done by transforming the order of samples in the sample data
    # frames to match the order of samples in the MOList object
    indexRNAMOList <- match(objMOList@RNAseqSamples$samples,
                             rownames(sampleDFRNAseq))
    indexSmallRNAMOList <- match(objMOList@SmallRNAseqSamples$samples, 
                                  rownames(sampleDFSmallRNAseq))
    sampleDFRNAseq <- sampleDFRNAseq[indexRNAMOList, ]
    sampleDFSmallRNAseq <- sampleDFSmallRNAseq[indexSmallRNAMOList, ]
    # Identify if the grouping variable is present in the sample data frames by
    # checking if the groupBy variable matches one column in the sample data
    # frames
    selRNA <- sapply(sampleDFRNAseq, function(x) all(x == groupByRNA))
    selSmallRNA <- sapply(sampleDFSmallRNAseq,
                          function(x) all(x == groupBySmallRNA))
    groupByVar <- intersect(names(selRNA)[selRNA], 
                            names(selSmallRNA)[selSmallRNA])
    if (length(groupByVar) == 0) {
      # Add the grouping variable to both sample data frames
      sampleDFRNAseq[, "groupBy"] <- groupByRNA
      sampleDFSmallRNAseq[, "groupBy"] <- groupBySmallRNA
    } else {
      # Rename the grouping variable to "groupBy" in both sample data frames
      colnames(sampleDFRNAseq)[colnames(sampleDFRNAseq) == 
                               groupByVar[1]] <- "groupBy"
      colnames(sampleDFSmallRNAseq)[colnames(sampleDFSmallRNAseq) ==
                                    groupByVar[1]] <- "groupBy"
    }

    objMOList$matchingRNAsmallRNA <- nnRNAMatch(sampleDFRNAseq, 
                                                      sampleDFSmallRNAseq)

  } else {
    # No additional sample information provided, use only the grouping variable
    # for matching
    sampleDFRNAseq <- data.frame(groupBy = groupByRNA)
    sampleDFSmallRNAseq <- data.frame(groupBy = groupBySmallRNA)
    objMOList$matchingRNAsmallRNA <- nnRNAMatch(sampleDFRNAseq, 
                                                      sampleDFSmallRNAseq)
  }
  return(objMOList)
}


# Use two ways to match: continuous grouping var to increase the weight of the grouping variable
# and binary grouping var to match separately between groups.












# Identify top genes that drive the differential expression of RNA or small RNAs
# top genes that drive the most variance separating the samples?

# use nearest neighbor matching to match the samples between the two groups
# then match the normalized counts of the top DE/variance genes between the two groups
# then use GENIE3 to predict the regulatory interactions between the top DE/variance genes and the small RNAs

# ensure when matching, the normalized counts of the small RNAs and the genes are generated from the same method
