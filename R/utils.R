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
#'
matchVecToMatrix <- function(vec, mat) {
  if (!is.null(vec) && !is.null(mat)) {
    if (length(vec) != ncol(mat)) {
      return(FALSE)
    } else {
      # Do nothing
    }
  } else {
    # Do nothing
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
  if (all(grepl(SRNA_SUFFIX, rownames(matchResult)))) {
    smallRNA <- rownames(matchResult)
    sampleMatch <- data.frame(
      smallRNA = gsub(SRNA_SUFFIX, "", smallRNA),
      RNA = matchResult[, 1]
    )
  } else {
    smallRNA <- matchResult[, 1]
    sampleMatch <- data.frame(
      smallRNA = gsub(SRNA_SUFFIX, "", smallRNA),
      RNA = rownames(matchResult)
    )
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
#' @return A data frame with the column with "colname" added
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


#' Format csAnno object to a data frame
#'
#' @keywords internal
#'
#' @param anno A csAnno object
#'
#' @return A data frame containing the annotation information
#'
#' @importFrom BiocGenerics as.data.frame
#'
csAnnoToDF <- function(anno) {
  annoDF <- BiocGenerics::as.data.frame(anno)
  colnames(annoDF)[1] <- "chr"
  annoDF <- annoDF %>%
    dplyr::select(
      chr, start, end, width, strand, Condition,
      annotation, ENSEMBL, SYMBOL, GENENAME
    ) %>%
    dplyr::arrange(chr, start)
  return(annoDF)
}


#' Get the names of enriched TFs
#'
#' @keywords internal
#'
#' @param enrichedMotifs A SummarizedExperiment object containing the motif
#'                       enrichment results
#' @return A vector of enriched TF names
#'
#' @importFrom dplyr %>%
#'
motifNames <- function(enrichedMotifs) {
  motifNames <- enrichedMotifs@elementMetadata@listData$motif.name

  # Some DNA binding elements are in the format of TF1::TF2 to indicate
  # dimer complexes. We extract both TF1 and TF2
  motifNames <- strsplit(motifNames, "::") %>%
    unlist() %>%
    unique()

  return(motifNames)
}


#' Extract a list of names respectively for down and up regulated genes
#'
#' @keywords internal
#'
#' @param deResult A data frame containing the differential expression results
#'
#' @return A list of names respectively for down and up regulated genes
#'
#' @importFrom dplyr %>% filter
#'
extractDirectionalGenes <- function(deResult) {
  downGenes <- deResult %>%
    dplyr::filter(logFC < 0) %>%
    rownames()
  upGenes <- deResult %>%
    dplyr::filter(logFC > 0) %>%
    rownames()
  return(list(down = downGenes, up = upGenes))
}


#' Annotate the type of gene
#'
#' @keywords internal
#'
#' @param geneName The name of a gene
#' @param annotation A list, with the name of each element defined to be the
#'                   type of gene, and the content of each element defined to be
#'                   a vector of gene names in that type
#'
#' @return The type of the gene
#'
findGeneType <- function(geneName, annotation) {
  geneType <- NULL
  for (i in seq_along(annotation)) {
    if (geneName %in% annotation[[i]]) {
      geneType <- names(annotation)[i]
      break
    } else {
      # Do nothing
    }
  }
  return(geneType)
}
findGeneType <- Vectorize(findGeneType, vectorize.args = "geneName")


# [END]
