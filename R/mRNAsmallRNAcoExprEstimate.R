# Purpose: Predictive estimation of regulatory small RNA - mRNA interaction
#          based on matched co-expression data
# Author: Jielin Yang
# Date: 2023-11-30
# Version: 1.0
# Bugs and Issues: None


# Global variable
SRNA_SUFFIX <- "_smallRNA"


#' Validate RNAseq and small RNAseq sample data frames for matching
#'
#' @keywords internal
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
    length(objMOList@smallRNAseqSamples$samples) !=
      nrow(sampleDFSmallRNAseq)) {
    stop("The number of samples in the sample data frames does not match the
      number of samples in the MOList object.")
  } else {
    # Continue
  }
  # Check for sample name matches
  if (!all(objMOList@RNAseqSamples$samples %in% rownames(sampleDFRNAseq)) ||
    !all(objMOList@smallRNAseqSamples$samples %in%
      rownames(sampleDFSmallRNAseq))) {
    stop("The sample names in the sample data frames do not match the sample
    names in the MOList object.")
  } else {
    # Continue
  }
  return(invisible(NULL))
}


#' Match RNAseq and small RNAseq samples based on continuous grouping variable
#'
#' @keywords internal
#'
#' @description This function performs optimal pair matching, a variant of
#'              nearest neighbor marching for RNAseq and small RNAseq samples
#'              based on a continuous grouping variable.
#'
#' @param sampleDF A data frame containing the sample information. Must follow
#'                 the following format:
#' \itemize{
#'  \item{row.names}{The sample names that correspond to the sample names in
#'                   the MOList object, with a suffix "_smallRNA" for small
#'                   RNAseq samples.}
#'  \item{groupBy}{The grouping variable for matching.}
#' }
#'
#' @return A data frame containing the matched sample information, with the
#'         following format:
#' \itemize{
#' \item{smallRNA}{The sample names of the small RNAseq samples. No suffix.}
#' \item{RNA}{The sample names of the RNAseq samples.}
#' }
#'
matchContinuous <- function(sampleDF) {
  # Use continuous matching algorithm
  matchVars <- colnames(sampleDF)[colnames(sampleDF) != "groupBy"]
  weight <- length(matchVars)

  # Define a "seq" variable to indicate whether the sample is an RNAseq sample
  # or a small RNAseq sample, which is also used for matching
  sampleDF <- labelSmallSizeGroup(sampleDF, SRNA_SUFFIX, "seq")

  # Performs weighted optimal pair matching
  if (length(matchVars) == 0) {
    # No additional variables for matching
    matchFormula <- seq ~ groupBy
  } else {
    # Additional variables for matching
    # Increase the weight of the grouping variable
    sampleDF[, paste0("groupBy", seq_len(weight))] <- sampleDF$groupBy
    matchVars <- colnames(sampleDF)[
      !(colnames(sampleDF) %in% c("groupBy", "seq"))
    ] # re-weighted
    matchFormula <- as.formula(paste0(
      "seq ~ groupBy + ",
      paste(matchVars, collapse = " + ")
    ))
  }
  # Perform matching using optimal pair
  matchResult <- MatchIt::matchit(
    matchFormula,
    data = sampleDF,
    method = "optimal",
    distance = "mahalanobis",
    ratio = 1
  )

  # Retrieve results and interpret match
  return(matchResultToDF(matchResult))
}


#' Match RNAseq and small RNAseq samples based on binary grouping variable
#'
#' @keywords internal
#'
#' @description This function performs optimal pair matching, a variant of
#'              nearest neighbor marching for RNAseq and small RNAseq samples
#'              based on a binary grouping variable. The matching is performed
#'              separately for the two groups.
#'
#' @param sampleDF A data frame containing the sample information.
#'                 Must follow the following format:
#' \itemize{
#' \item{row.names}{The sample names that correspond to the sample names in
#'                  the MOList object, with a suffix "_smallRNA" for small
#'                  RNAseq samples.}
#' \item{groupBy}{The grouping variable for matching, must be binary.}
#' }
#'
#' @return A data frame containing the matched sample information, with the
#'        following format:
#' \itemize{
#' \item{smallRNA}{The sample names of the small RNAseq samples. No suffix.}
#' \item{RNA}{The sample names of the RNAseq samples.}
#' }
#'
#' @importFrom MatchIt matchit
#' @import optmatch
#'
matchBinary <- function(sampleDF) {
  # Covariates for matching
  matchVars <- colnames(sampleDF)[colnames(sampleDF) != "groupBy"]

  # Matching formula
  matchFormula <- as.formula(paste0(
    "seq ~ ",
    paste(matchVars, collapse = " + ")
  ))

  # Separate the sample data frame into two data frames based on the grouping
  # variable
  group1 <- sampleDF$groupBy == unique(sampleDF$groupBy)[1]
  sampleDF1 <- sampleDF[group1, ] %>% labelSmallSizeGroup(SRNA_SUFFIX, "seq")
  sampleDF2 <- sampleDF[!group1, ] %>% labelSmallSizeGroup(SRNA_SUFFIX, "seq")

  # Perform matching using optimal pair
  matchResult1 <- MatchIt::matchit(
    matchFormula,
    data = sampleDF1,
    method = "optimal",
    distance = "mahalanobis",
    ratio = 1
  )
  matchResult2 <- MatchIt::matchit(
    matchFormula,
    data = sampleDF2,
    method = "optimal",
    distance = "mahalanobis",
    ratio = 1
  )
  matchedSampleDF <- rbind(
    matchResultToDF(matchResult1),
    matchResultToDF(matchResult2)
  )
  return(matchedSampleDF)
}


#' Performing optimal pair matching for RNAseq and small RNAseq samples
#'
#' @description This function performs optimal pair matching, a variant of
#'              nearest neighbour marching for RNAseq and
#'              small RNAseq samples using two sample data frames. The two data
#'              frames should contain the same number of variables, including
#'              the main grouping variable for DE analysis.
#'
#' @details Depending on the nature of the grouping variables, either continuous
#'          or binary, the matching algorithm will be different. For continuous
#'          grouping variables, the grouping variable will has a weight equal to
#'          the number of other variables in the matching process. This is to
#'          ensure that the algorithm matches samples with the same grouping
#'          variable first, and then matches samples with similar values of the
#'          additional parameters. For binary grouping variables, the algorithm
#'          forces the matching process to be within the same group, and then
#'          performs optimal pair matching on other parameters.
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
#' @return A list of two numeric vectors with equal length, where each element
#'         in the first vector is the index of the RNAseq sample, and the
#'         corresponding element in the second vector is the index of the
#'         matched small RNAseq sample.
#'
#' @importFrom MatchIt matchit
#' @import optmatch
#'
nnRNAMatch <- function(objMOList, sampleDFRNAseq, sampleDFSmallRNAseq) {
  # Check if the two sample data frames have the same number of variables
  if (ncol(sampleDFRNAseq) != ncol(sampleDFSmallRNAseq)) {
    stop("The two sample data frames do not have the same number of variables.")
  } else {
    # Continue
  }

  # Concatenate the two sample data frames
  rownames(sampleDFSmallRNAseq) <- paste0(
    rownames(sampleDFSmallRNAseq),
    SRNA_SUFFIX
  )
  sampleDF <- rbind(sampleDFRNAseq, sampleDFSmallRNAseq)

  if (length(unique(sampleDFRNAseq$groupBy)) > 2) {
    sampleMatchDF <- matchContinuous(sampleDF)
  } else {
    sampleMatchDF <- matchBinary(sampleDF)
  }

  # Convert the sample names to indices
  indexRNAseq <- match(
    sampleMatchDF$RNA,
    objMOList@RNAseqSamples$samples
  )
  indexSmallRNAseq <- match(
    sampleMatchDF$smallRNA,
    objMOList@smallRNAseqSamples$samples
  )
  return(list(indexRNAseq = indexRNAseq, indexSmallRNAseq = indexSmallRNAseq))
}



#' Match RNAseq and small RNAseq samples
#'
#' @description This function generates a one-to-one matching between the
#'              RNAseq samples and the small RNAseq samples using the
#'              optimal pair matching algorithm, which is a variant of nearest
#'              neighbor matching. The matching is based on the grouping
#'              variable and other variables in the sample data frames.
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
#'         in the RNAseqSamples and smallRNAseqSamples slots of the MOList
#'         object. The slot is a list under the following structure:
#' \itemize{
#' \item{indexRNAseq}{A numeric vector containing the indices of the RNAseq
#'                    samples.}
#' \item{indexSmallRNAseq}{A numeric vector containing the indices of the small
#'                         RNAseq samples.}
#' }
#'
#' @importFrom MatchIt matchit
#' @import optmatch
#'
#' @export
#'
#' @examples
#' # Suppose we have an MOList object with RNAseq and small RNAseq data
#'
#' # Example 1: No additional sample information provided
#' \dontrun{
#' # Perform matching
#' myMOList <- matchSamplesRNAsmallRNA(myMOList)
#' }
#'
#' # Example 2: Additional sample information provided
#' \dontrun{
#' # Suppose we have two data frames containing the sample information
#' sampleDFRNAseq <- data.frame(
#'   Age = c(1, 2, 3, 4, 5),
#'   Sex = c("M", "M", "F", "F", "F")
#' )
#' sampleDFSmallRNAseq <- data.frame(
#'   Age = c(1, 2, 3, 4, 5),
#'   Sex = c("M", "M", "F", "F", "F")
#' )
#' # Perform matching
#' myMOList <- matchSamplesRNAsmallRNA(myMOList,
#'   sampleDFRNAseq = sampleDFRNAseq,
#'   sampleDFSmallRNAseq = sampleDFSmallRNAseq
#' )
#' # Will account for both age and sex
#' }
#'
#' # Example 3: Additional sample information provided, with matching variables
#' \dontrun{
#' # Suppose we just want to match based on age
#' myMOList <- matchSamplesRNAsmallRNA(myMOList,
#'   sampleDFRNAseq = sampleDFRNAseq,
#'   sampleDFSmallRNAseq = sampleDFSmallRNAseq,
#'   varMatch = "Age"
#' )
#' # Will only account for age
#' }
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
  groupBySmallRNA <- objMOList@smallRNAseqSamples$groupBy

  # Only use sample information when both sampleDFRNAseq and
  # sampleDFSmallRNAseq are provided
  if (!is.null(sampleDFRNAseq) && !is.null(sampleDFSmallRNAseq)) {
    validateSampleDFs(objMOList, sampleDFRNAseq, sampleDFSmallRNAseq)
    if (!is.null(varMatch)) {
      invalidVars <- setdiff(varMatch, intersect(
        colnames(sampleDFRNAseq),
        colnames(sampleDFSmallRNAseq)
      ))
      if (length(invalidVars) > 0) {
        warning(paste("The following variables are not present in both sample
                      data frames and are ignored:", invalidVars))
        varMatch <- setdiff(varMatch, invalidVars)
      } else {
        # Continue
      }
    } else {
      varMatch <- intersect(
        colnames(sampleDFRNAseq),
        colnames(sampleDFSmallRNAseq)
      )
    }
    # Keep only the usable sample information
    sampleDFRNAseq <- sampleDFRNAseq[, varMatch]
    sampleDFSmallRNAseq <- sampleDFSmallRNAseq[, varMatch]

    # Since the key grouping variable may or may not be present in the selected
    # matching variables, we need to check if it is present, and if not, add it
    # to the matching variables (and the sample data frames)
    # This is done by transforming the order of samples in the sample data
    # frames to match the order of samples in the MOList object
    indexRNAMOList <- match(
      objMOList@RNAseqSamples$samples,
      rownames(sampleDFRNAseq)
    )
    indexSmallRNAMOList <- match(
      objMOList@smallRNAseqSamples$samples,
      rownames(sampleDFSmallRNAseq)
    )
    sampleDFRNAseq <- sampleDFRNAseq[indexRNAMOList, ]
    sampleDFSmallRNAseq <- sampleDFSmallRNAseq[indexSmallRNAMOList, ]
    # Identify if the grouping variable is present in the sample data frames by
    # checking if the groupBy variable matches one column in the sample data
    # frames
    selRNA <- sapply(sampleDFRNAseq, function(x) all(x == groupByRNA))
    selSmallRNA <- sapply(
      sampleDFSmallRNAseq,
      function(x) all(x == groupBySmallRNA)
    )
    groupByVar <- intersect(
      names(selRNA)[selRNA],
      names(selSmallRNA)[selSmallRNA]
    )
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

    objMOList$matchingRNAsmallRNA <- nnRNAMatch(
      objMOList,
      sampleDFRNAseq,
      sampleDFSmallRNAseq
    )
  } else {
    # No additional sample information provided, use only the grouping variable
    # for matching
    sampleDFRNAseq <- data.frame(groupBy = groupByRNA)
    sampleDFSmallRNAseq <- data.frame(groupBy = groupBySmallRNA)
    objMOList$matchingRNAsmallRNA <- nnRNAMatch(
      objMOList,
      sampleDFRNAseq,
      sampleDFSmallRNAseq
    )
  }
  return(objMOList)
}


#' Export the match result of RNAseq and small RNAseq samples
#'
#' @description This function exports the match result of RNAseq and small
#'              RNAseq samples.
#'
#' @param objMOList An object of class MOList
#'
#' @return A data frame containing the match result, with the following format:
#' \itemize{
#' \item{RNAseq}{The sample names of the RNAseq samples.}
#' \item{SmallRNAseq}{The sample names of the small RNAseq samples.}
#' }
#'
#' @export
#'
#' @examples
#' # Suppose we have an MOList object with RNAseq and small RNAseq data
#' \dontrun{
#' # Perform matching
#' myMOList <- matchSamplesRNAsmallRNA(myMOList)
#' # Export the match result
#' exportMatchResult(myMOList)
#' }
#'
exportMatchResult <- function(objMOList) {
  if (is.null(objMOList$matchingRNAsmallRNA)) {
    stop("The MOList object does not contain matching information.")
  } else {
    # Continue
  }
  matchResult <- data.frame(
    RNAseq = objMOList@RNAseqSamples$samples[
      objMOList$matchingRNAsmallRNA$indexRNAseq
    ],
    SmallRNAseq = objMOList@smallRNAseqSamples$samples[
      objMOList$matchingRNAsmallRNA$indexSmallRNAseq
    ]
  )
  return(matchResult)
}


#' Run GENIE3 for predicted interactions
#'
#' @keywords internal
#'
#' @description This function runs GENIE3 for predicted interactions.
#'
#' @param exprMatrix A numeric matrix containing the expression data.
#' @param regulators A character vector containing the names of the regulators.
#' @param ntree The number of trees to be used in the ensemble learning.
#' @param nthreads The number of threads to be used in the ensemble learning.
#' @param treeMethod The tree method to be used in the ensemble learning, either
#'                   "RF" for random forest or "ET" for extra trees. See the
#'                   documentation of the GENIE3 package for details.
#' @param seed The seed to be used in the ensemble learning for reproducibility.
#'
#' @return A list containing the predicted interactions, with the
#'        following format:
#' \itemize{
#' \item{regulator}{The names of the small RNAs that are predicted to regulate
#'                       the target genes.}
#' \item{target}{The names of the target genes that are predicted to be
#'                   regulated by the small RNAs.}
#' \item{weight}{The weight of the predicted interactions.}
#' }
#'
#' @importFrom GENIE3 GENIE3 getLinkList
#'
runGENIE3 <- function(exprMatrix,
                      regulators,
                      ntree = 1000,
                      nthreads = 1,
                      treeMethod = "RF",
                      seed = 91711) {
  cat("Running GENIE3 for predicted interactions...\n")
  cat("  Number of trees: ", ntree, "\n")
  cat("  Number of threads: ", nthreads, "\n")
  cat("  Tree method: ", treeMethod, "\n")
  cat("This may take a while...\n")

  # Set seed for reproducibility
  set.seed(seed)

  # Run GENIE3
  weightedMatrix <- GENIE3::GENIE3(
    exprMatrix,
    regulators = regulators,
    nTrees = ntree,
    nCores = nthreads,
    treeMethod = treeMethod
  )
  set.seed(NULL)

  # The weighted matrix contains a weight for each regulator-target pair
  weightedAdjList <- GENIE3::getLinkList(weightedMatrix)

  # Remove any predicted target genes which are also regulators
  weightedAdjList <- weightedAdjList %>% dplyr::filter(
    !(targetGene %in% regulators)
  )

  # Consider only the top 20% of the predicted interactions
  weightedAdjList <- weightedAdjList %>% dplyr::filter(
    weight >= quantile(weight, 0.8)
  )

  # Convert the weighted adjacency list to a list of three vectors
  weightedAdjList <- list(
    regulator = weightedAdjList[, 1],
    target = weightedAdjList[, 2]
  )

  rm(weightedMatrix)

  return(weightedAdjList)
}


#' Predict small RNA - mRNA interactions via tree-based ensemble learning
#'
#' @keywords internal
#'
#' @description This function predicts small NC RNA - mRNA interactions via
#'              tree-based ensemble learning using the GENIE3 algorithm.
#'              The prediction is directed from small RNA to mRNA.
#'
#' @param mRNATopTag A TOPTag object containing the top mRNAs that are
#'                   differentially expressed cross the conditions.
#' @param smallRNATopTag A TOPTag object containing the top small RNAs that are
#'                       differentially expressed cross the conditions.
#' @param smallRNATypes A character vector containing the types of small RNAs to
#'                      be used for prediction. The types should be a valid type
#'                      defined during the annotation process. Must be a subset
#'                      of the types: "miRNA", "piRNA", "tRNA", "snRNA",
#'                      "snoRNA", "circRNA", or a single string "all" to use all
#'                      types.
#' @param annoSncRNA A list containing the transcripts of the small RNAs
#'                   in each respective type. The list should be named
#'                   consistently as defined in annotateSmallRNA().
#' @param matchingRNAsmallRNA A list of two numeric vectors with equal length,
#'                            where each element in the first vector is the
#'                            index of the RNAseq sample, and the corresponding
#'                            element in the second vector is the index of the
#'                            matched small RNAseq sample. These indices match
#'                            to the sample names in the RNAseqSamples and
#'                            smallRNAseqSamples slots of the MOList object.
#' @param ntree The number of trees to be used in the ensemble learning.
#' @param nthreads The number of threads to be used in the ensemble learning.
#' @param treeMethod The tree method to be used in the ensemble learning, either
#'                   "RF" for random forest or "ET" for extra trees. See the
#'                   documentation of the GENIE3 package for details.
#' @param seed The seed to be used in the ensemble learning for reproducibility.
#'
#' @return A list containing the predicted interactions, with the
#'         following format:
#' \itemize{
#' \item{regulator}{The names of the small RNAs that are predicted to regulate
#'                  the target genes.}
#' \item{target}{The names of the target genes that are predicted to be
#'               regulated by the small RNAs.}
#' \item{weight}{The weight of the predicted interactions.}
#' }
#'
#' @note The final result includes all predicted interactions, without any
#'       filtering on the weight of the predicted interactions.
#'
#' @importFrom GENIE3 GENIE3 getLinkList
#'
#'
predictSmallRNAmRNAcoExpr <- function(mRNATopTag,
                                      smallRNATopTag,
                                      smallRNATypes = c(
                                        "miRNA", "piRNA",
                                        "snRNA", "snoRNA",
                                        "circRNA", "tRNA"
                                      ),
                                      annoSncRNA,
                                      matchingRNAsmallRNA,
                                      ntree = 1000,
                                      nthreads = 1,
                                      treeMethod = "RF",
                                      seed = 91711) {
  # Check if required annotation and matching information are provided
  if (is.null(annoSncRNA) || length(annoSncRNA) == 0) {
    stop("The annotation information of small RNAs is not provided.")
  } else if (is.null(matchingRNAsmallRNA)) {
    stop("The matching information of RNAseq and small RNAseq samples is not
      provided.")
  } else {
    # Do nothing
  }

  # Show the user about construction information
  cat("Constructing a combined expression matrix for GENIE3...\n")
  cat("Criteria for defining the top differentially expressed genes:\n")
  cat("  log2 fold change cutoff: ", mRNATopTag@logFCCutoff, "\n")
  cat("  adjusted p-value cutoff: ", mRNATopTag@pCutoff, "\n")
  if (mRNATopTag@topGenes > 1) {
    cat("  top", mRNATopTag@topGenes, " genes selected\n")
  } else {
    cat("  top", mRNATopTag@topGenes * 100, "% genes selected\n")
  }
  cat("Criteria for defining the top differentially expressed small RNAs:\n")
  cat("  log2 fold change cutoff: ", smallRNATopTag@logFCCutoff, "\n")
  cat("  adjusted p-value cutoff: ", smallRNATopTag@pCutoff, "\n")
  if (smallRNATopTag@topGenes > 1) {
    cat("  top", smallRNATopTag@topGenes, " small RNAs selected\n")
  } else {
    cat("  top", smallRNATopTag@topGenes * 100, "% small RNAs selected\n")
  }

  # Construct a master expression matrix
  exprMatrixRNA <- mRNATopTag@normalizedCounts %>%
    as.data.frame() %>%
    dplyr::select(matchingRNAsmallRNA$indexRNAseq)
  exprMatrixSmallRNA <- smallRNATopTag@normalizedCounts %>%
    as.data.frame() %>%
    dplyr::select(matchingRNAsmallRNA$indexSmallRNAseq)
  sampleNames <- paste0("match_", seq_len(ncol(exprMatrixRNA)))
  colnames(exprMatrixRNA) <- sampleNames
  colnames(exprMatrixSmallRNA) <- sampleNames
  exprMatrix <- rbind(exprMatrixRNA, exprMatrixSmallRNA) %>% as.matrix()

  # Obtain the set of regulators for GENIE3
  regulators <- unlist(annoSncRNA[smallRNATypes]) %>%
    intersect(rownames(exprMatrixSmallRNA))

  # Run GENIE3 for predicted interactions
  weightedAdjList <- runGENIE3(
    exprMatrix,
    regulators = regulators,
    ntree = ntree,
    nthreads = nthreads,
    treeMethod = treeMethod,
    seed = seed
  )

  return(weightedAdjList)
}
