# Purpose: Performing differential analysis on the omics data
# Author: Jielin Yang
# Date: 2023-10-30
# Version: 1.0
# Bugs and Issues: None


# Define a global variable for the count-based omics data
RNA <- "RNAseq"
SMALLRNA <- "smallRNAseq"
PROTEIN <- "proteomics"
ATAC <- "ATACseq"
COUNT_OMICS <- c(RNA, SMALLRNA, PROTEIN)


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
validateDataAnno <- function(objMOList, annoList) {
  # Validate the correctness of the MOList object
  validateMOList(objMOList)
  # Check of annotations on the samples
  if (!matchVecToMatrix(annoList$RNAseq, objMOList@RNAseq) ||
    !matchVecToMatrix(annoList$smallRNAseq, objMOList@smallRNAseq) ||
    !matchVecToMatrix(annoList$proteomics, objMOList@proteomics)) {
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
#' \itemize{
#' \item code(RNAseq): The RNAseq data in the MOList object is filtered and
#'                     updated, if omic = "RNAseq"
#' \item code(smallRNAseq): The small RNAseq data in the MOList object is
#'                           filtered and updated, if omic = "smallRNAseq"
#' \item code(proteomics): The protein data in the MOList object is filtered and
#'                         updated, if omic = "proteomics"
#' \item code{ATACpeaks}: A list containing the ATAC peaks for condition 1 and
#'                        condition 2, kept unchanged
#' }
#'
#' @importFrom edgeR cpm
#'
#' @references
#' \insertRef{robinson2010edger}{IntegraTRN}
#'
filterGeneCounts <- function(objMOList, omic) {
  if (!(omic %in% COUNT_OMICS)) {
    stop("The input omics data is not supported for filtering")
  }

  omicData <- getRawData(objMOList, omic)
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
      proteomics = proteomics(objMOList) <- omicData
    )
    return(objMOList)
  }
}


#' Perform differential expression analysis on the count-based omics data using
#' DESeq2
#'
#' @keywords internal
#'
#' @description This function performs differential expression analysis on the
#'              count-based omics data using DESeq2. The RNAseq, small RNAseq,
#'              and protein data are supported for differential expression
#'              analysis.
#'
#' @param filteredCounts A numeric matrix containing the filtered count-based
#'                     omics data
#' @param groupBy A vector specifying the grouping information for the omics
#'                data, must be the same length as the number of samples in the
#'                omics data
#' @param batch A character vector specifying the batch information for the
#'              omics data, must be the same length as the number of samples in
#'              the omics data, used for batch correction. Can be NULL if no
#'              batch correction is needed.
#'
#' @return A DETag object containing the differential expression analysis
#'         results, and the method DESeq2
#' \itemize{
#' \item \code{DEResult}: A data frame containing the differential expression
#'                        analysis results
#' \item \code{method}: The character string "DESeq2"
#' }
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results counts
#' @importFrom dplyr %>%
#'
#' @references
#' \insertRef{love2014moderated}{IntegraTRN}
#'
diffExprDESeq2 <- function(filteredCounts, groupBy, batch = NULL) {
  # Generate DESeqDataSet as the base object
  designList <- DESeqDesign(groupBy, batch)
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = filteredCounts,
    colData = designList$colData,
    design = designList$design
  )

  # Perform differential expression analysis and obtain results
  dds <- DESeq2::DESeq(dds)
  DEResult <- DESeq2::results(dds) %>% as.data.frame()
  normalizedCounts <- DESeq2::counts(dds, normalized = TRUE) %>% as.matrix()

  # Create the DETag object
  deTag <- DETag(
    DEResult = DEResult,
    method = DESEQ2,
    normalizedCounts = normalizedCounts
  )

  # Free up memory for large operational objects
  rm(dds)

  # Return the results
  return(deTag)
}


#' Perform differential expression analysis on the count-based omics data using
#' EdgeR
#'
#' @keywords internal
#'
#' @description This function performs differential expression analysis on the
#'              count-based omics data using EdgeR. The RNAseq, small RNAseq,
#'              and protein data are supported for differential expression
#'              analysis.
#'
#' @param filteredCounts A numeric matrix containing the filtered count-based
#'                     omics data
#' @param groupBy A vector specifying the grouping information for the omics
#'                data, must be the same length as the number of samples in the
#'                omics data
#' @param batch A character vector specifying the batch information for the
#'              omics data, must be the same length as the number of samples in
#'              the omics data, used for batch correction. Can be NULL if no
#'              batch correction is needed.
#'
#' @return A DETag object containing the differential expression analysis
#'         results, and the method DESeq2
#' \itemize{
#' \item \code{DEResult}: A data frame containing the differential expression
#'                        analysis results
#' \item \code{method}: The character string "EdgeR"
#' }
#'
#' @importFrom edgeR DGEList estimateDisp calcNormFactors glmQLFit glmQLFTest
#' @importFrom edgeR topTags cpm
#' @importFrom dplyr %>%
#'
#' @references
#' \insertRef{robinson2010edger}{IntegraTRN}
#'
diffExprEdgeR <- function(filteredCounts, groupBy, batch = NULL) {
  # Generate DGEList as the base object
  dge <- edgeR::DGEList(counts = filteredCounts, group = groupBy)

  # Construct GLM model
  if (!is.null(batch)) {
    design <- stats::model.matrix(~ batch + groupBy)
  } else {
    design <- stats::model.matrix(~groupBy)
  }

  # Normalization and dispersion estimation
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateDisp(dge, design = design)

  # Perform differential expression analysis and obtain results
  fit <- edgeR::glmQLFit(dge, design = design)
  qlf <- edgeR::glmQLFTest(fit, coef = ncol(design))
  DEResult <- edgeR::topTags(qlf, n = nrow(filteredCounts)) %>% as.data.frame()
  normalizedCounts <- edgeR::cpm(dge) %>% as.matrix()

  # Create the DETag object
  deTag <- DETag(
    DEResult = DEResult,
    method = EDGER,
    normalizedCounts = normalizedCounts
  )

  # Free up memory for large operational objects
  rm(dge, fit, qlf)

  # Return the results
  return(deTag)
}


#' Differential expression analysis of count-based omics data
#'
#' @keywords internal
#'
#' @description This function performs differential expression analysis on the
#'              count-based omics data. The RNAseq, small RNAseq, and protein
#'              data are supported for differential expression analysis.
#'              Calling this function again on the same omics data will
#'              overwrite the previous results.
#'
#' @note Preconditions:
#' \itemize{
#' \item The input omics data must be a MOList object
#' \item The input omics data must be filtered by filterGeneCounts
#' \item The batch information must have the same length as the number of
#'       samples in the omics data
#' }
#'
#' @param objMOList A MOList object containing the omics data
#' @param omic A character string specifying the omics data to be analyzed
#'             must be one of "RNAseq", "smallRNAseq", and "proteomics"
#' @param batch A character vector specifying the batch information for the
#'              omics data, must be the same length as the number of samples in
#'              the omics data, used for batch correction. Can be NULL if no
#'              batch correction is needed.
#' @param program A character string specifying the program used for the
#'               analysis, DESeq2 or EdgeR
#'
#' @return An MOList object containing the differential expression analysis
#'         results. Results are appended to the original MOList object as
#'         list elements named "DERNAseq", "DEsmallRNAseq", and "DEProteomics"
#'         for RNAseq, small RNAseq, and protein data, respectively.
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
#' \item \code{DERNAseq}: A DETag object containing the differential
#'                        expression analysis results for the RNAseq data
#' \item \code{DEsmallRNAseq}: A DETag object containing the differential
#'                             expression analysis results for the smallRNAseq
#'                             data
#' \item \code{DEproteomics}: A DETag object containing the differential
#'                            expression analysis results for the proteomics
#'                            data
#' }
#'
#' @importFrom dplyr %>%
#' @importFrom edgeR DGEList estimateDisp calcNormFactors glmQLFit glmQLFTest
#'             topTags
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#'
#' @references
#' \insertRef{love2014moderated}{IntegraTRN}
#'
#' \insertRef{robinson2010edger}{IntegraTRN}
#'
countDiffExpr <- function(objMOList, omic, batch, program = DESEQ2) {
  if (!(omic %in% COUNT_OMICS)) {
    stop("The input omics data is not supported for differential analysis")
  } else if (is.null(getRawData(objMOList, omic))) {
    # Nothing to do, return the original object
    return(objMOList)
  } else {
    # Do nothing
  }

  cat("Performing differential analysis on", omic, "data\n\n")

  # Differential analysis based on selected program
  filtedCounts <- getRawData(objMOList, omic)
  if (program == DESEQ2) {
    DEResult <- diffExprDESeq2(
      filteredCounts = filtedCounts,
      groupBy = getSampleInfo(objMOList, omic)$groupBy,
      batch = batch
    )
  } else if (program == EDGER) {
    DEResult <- diffExprEdgeR(
      filteredCounts = filtedCounts,
      groupBy = getSampleInfo(objMOList, omic)$groupBy,
      batch = batch
    )
  } else {
    stop("The input program is not supported for differential analysis")
  }

  # Update the MOList object
  objMOList[[paste0("DE", omic)]] <- DEResult
  return(objMOList)
}


#' Perform differential analysis on the omics data
#'
#' @aliases diffOmics
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
#' @param program A character string specifying the program used for the
#'                analysis, DESeq2 or EdgeR
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom edgeR DGEList estimateDisp calcNormFactors glmQLFit glmQLFTest
#'             topTags
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom GenomicRanges makeGRangesFromDataFrame reduce
#'
#' @return An MOList object containing the differential analysis results
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
#' \item \code{DERNAseq}: A DETag object containing the differential
#'                        expression analysis results for the RNAseq data
#' \item \code{DEsmallRNAseq}: A DETag object containing the differential
#'                             expression analysis results for the smallRNAseq
#'                             data
#' \item \code{DEproteomics}: A DETag object containing the differential
#'                            expression analysis results for the proteomics
#'                            data
#' \item \code{DEATAC}: A DETag object object containing the differential
#'                      accessible accessible regions for the ATACseq data
#' }
#'
#' @references
#' \insertRef{love2014moderated}{IntegraTRN}
#'
#' \insertRef{robinson2010edger}{IntegraTRN}
#'
#' \insertRef{lawrence2013software}{IntegraTRN}
#'
#' @examples
#' # Generate example datasets
#' rnaseq <- matrix(sample(0:100, 1000, replace = TRUE), nrow = 100, ncol = 10)
#' rnaseqSGroupBy <- rep(c("A", "B"), each = 5)
#' smallRNAseq <- matrix(sample(0:100, 1000, replace = TRUE),
#'   nrow = 100, ncol = 5
#' )
#' smallRNAseqSGroupBy <- rep(c("A", "B"), each = 3)[1:5]
#' protein <- matrix(sample(0:100, 1000, replace = TRUE), nrow = 100, ncol = 5)
#' proteinSGroupBy <- rep(c("A", "B"), each = 3)[1:5]
#'
#' # Create a MOList object
#' objMOList <- MOList(
#'   RNAseq = rnaseq, RNAGroupBy = rnaseqSGroupBy,
#'   smallRNAseq = smallRNAseq, smallRNAGroupBy = smallRNAseqSGroupBy,
#'   proteomics = protein, proteomicsGroupBy = proteinSGroupBy
#' )
#'
#' # Perform differential analysis
#' objMOList <- diffOmics(objMOList)
#'
diffOmics <- function(objMOList,
                      rnaseqBatch = NULL,
                      smallRnaBatch = NULL,
                      proteinBatch = NULL,
                      program = "DESeq2") {
  # Validating inputs
  validateDataAnno(objMOList, list(
    RNAseq = rnaseqBatch,
    smallRNAseq = smallRnaBatch,
    proteomics = proteinBatch
  ))

  # Data processing and differential expression of count-based omics data
  # Here filterGeneCounts and diffExprCount are called on all types of
  # count-based omics data, even if they could be NULL
  for (omic in COUNT_OMICS) {
    # Pass if the omics data is NULL
    if (length(getRawData(objMOList, omic)) == 0) {
      # Do nothing
    } else {
      # Filter raw counts
      objMOList <- filterGeneCounts(objMOList, omic)
      # Differential expression
      objMOList <- countDiffExpr(
        objMOList = objMOList,
        omic = omic,
        batch = switch(omic,
          RNAseq = rnaseqBatch,
          smallRNAseq = smallRnaBatch,
          proteomics = proteinBatch
        ),
        program = program
      )
    }
  }

  # Construct a master list of differentially accessible regions with merged
  # intra-condition peaks
  objMOList <- processPeakOverlap(objMOList)

  return(objMOList)
}


# [END]
