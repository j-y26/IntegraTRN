# Purpose: Performing differential analysis on the omics data
# Author: Jielin Yang
# Date: 2023-10-30
# Version: 1.0
# Bugs and Issues: Currently, only DESeq2 is supported for differential analysis


# Define global variables
DESEQ2 <- "DESeq2"
EDGER <- "edgeR"
ATAC_GRANGE <- "GRangesATAC"
METHODS <- c(DESEQ2, EDGER, ATAC_GRANGE)
COUNT_DEFIELDS <- c("logFC", "pvalue", "padj")
EDGER_FIELDS <- c("logFC", "PValue", "FDR")
DESEQ2_FIELDS <- c("log2FoldChange", "pvalue", "padj")


#' @name DETag-class
#' @title DETag S4 class
#' @aliases DETag-class
#'
#' @description This class is used to store the differential expression analysis
#'              results for the count-based omics data. Although the
#'              functionality of the class itself can be easily implemented
#'              by using other data structures, here use specifically designed
#'              S4 class to allow future extension of the program, e.g.
#'              providing different methods for differential expression
#'              analysis while maintaining a consistent interface.
#'
#' @slot DEResult A data frame containing the differential expression analysis
#'                results, with each row representing a gene and each column
#'                representing a feature. The exact columns depend on the
#'                method used for the analysis
#' @slot method A character string specifying the program used for the analysis,
#'              must be one of the valid METHODS defined
#' @slot normalizedCounts A matrix containing the normalized counts for the
#'                        genes, with each row representing a gene and each
#'                        column representing a sample
#'
#' @exportClass DETag
#'
#' @references
#' \insertRef{robinson2010edger}{IntegraTRN}
#'
#' \insertRef{love2014moderated}{IntegraTRN}
#'
setClass("DETag", slots = c(
  DEResult = "data.frame",
  method = "character",
  normalizedCounts = "matrix"
))


#' Validator of the DETag class slots
#'
#' @keywords internal
#'
#' @description This function validates the DETag class
#'
#' @param DEResult A data frame containing the differential expression analysis
#'                 results
#' @param method A character string specifying the program used for the
#'               analysis, must be one of the valid METHODS defined
#'
#' @return NULL
#'
validateDETagSlots <- function(DEResult, method, normalizedCounts = NULL) {
  # Validate each slot
  # DEResult
  if (!inherits(DEResult, "data.frame")) {
    stop("The DEResult slot must be a data frame")
  } else if (any(is.na(DEResult[, 1:3]))) {
    stop("Some key differential analysis results are missing")
  } else {
    # Do nothing
  }
  # method
  if (!is.character(method) || !method %in% METHODS) {
    stop("The method slot must be a valid character string, see ?DETag for
         more information")
  } else {
    # Do nothing
  }
  # normalizedCounts
  if (!is.null(normalizedCounts) &&
    (!is.matrix(normalizedCounts) || !is.numeric(normalizedCounts))) {
    stop("The normalizedCounts slot must be a numeric matrix")
  } else {
    # Do nothing
  }

  return(invisible(NULL))
}


#' Constructor for the DETag class
#'
#' @aliases DETag
#'
#' @description This function creates an object of the DETag class
#'
#' @param DEResult A data frame containing the differential expression analysis
#'                 results
#' @param method A character string specifying the program used for the analysis
#'
#' @param normalizedCounts A matrix containing the normalized counts for the
#'                         genes, with each row representing a gene and each
#'                         column representing a sample
#'
#' @return An object of the DETag class
#' \itemize{
#' \item \code{DEResult}: A data frame containing the differential expression
#'                        analysis results
#' \item \code{method}: A character string specifying the program used for the
#'                      analysis
#' }
#'
#' @export
#'
#' @examples
#' # Create an example data frame
#' deResult <- data.frame(
#'   gene = paste0("gene", seq_len(10)),
#'   logFC = runif(10),
#'   adj.P.Val = runif(10)
#' )
#'
#' # Create a normalized counts matrix
#' normalizedCounts <- matrix(runif(100), nrow = 10)
#'
#' # Create an object of the DETag class
#' deTag <- DETag(deResult, "DESeq2", normalizedCounts)
#'
#' # Check the class of the object
#' class(deTag)
#'
DETag <- function(DEResult, method, normalizedCounts = NULL) {
  # Validate the input
  validateDETagSlots(DEResult, method, normalizedCounts)
  # Create the object
  if (method %in% c(DESEQ2, EDGER)) {
    # Count-based differential analysis
    newDETag <- new("DETag",
      DEResult = DEResult,
      method = method,
      normalizedCounts = normalizedCounts
    )
  } else if (method == ATAC_GRANGE) {
    # Differential accessibility analysis
    newDETag <- new("DETag", DEResult = DEResult, method = method)
  } else {
    stop("The method is not supported")
  }
  # Return the object
  return(newDETag)
}


# Define the S4 method for the DETag class

#' Print DETag
#'
#' @aliases print,DETag-method
#'
#' @description This function prints the DETag object
#'
#' @param object An object of the DETag class
#'
#' @return NULL
#'
#' @export
#'
#' @examples
#' # Create an example data frame
#' deResult <- data.frame(
#'   gene = paste0("gene", seq_len(10)),
#'   logFC = runif(10),
#'   adj.P.Val = runif(10)
#' )
#'
#' # Create a normalized counts matrix
#' normalizedCounts <- matrix(runif(100), nrow = 10)
#'
#' # Create an object of the DETag class
#' deTag <- DETag(deResult, "DESeq2", normalizedCounts)
#'
#' # Print the object
#' print(deTag)
#'
#' # or simply do
#' deTag
#'
setMethod("show", "DETag", function(object) {
  cat("A DETag S4 object\n")
  cat("Method for differential analysis: ", object@method, "\n\n")
  # Some more DE details for count-based DE
  if (object@method == DESEQ2) {
    cat("Number of genes: ", nrow(object@DEResult), "\n")
    cat(
      "Number of genes with adjusted p-value < 0.05: ",
      sum(object@DEResult$padj < 0.05, na.rm = TRUE), "\n\n"
    )
  } else if (object@method == EDGER) {
    cat("Number of genes: ", nrow(object@DEResult), "\n")
    cat(
      "Number of genes with FDR < 0.05: ",
      sum(object@DEResult$FDR < 0.05, na.rm = TRUE), "\n\n"
    )
  } else {
    # Do nothing
  }
  cat("A snap shot for differential analysis results:\n")
  print(utils::head(object@DEResult))
  cat("\n")
  cat("To access the full results, use the exportDE S4 method\n")
  invisible(object)
})


#' Export the differential analysis results from a DETag object
#'
#' @aliases exportDE,DETag-method
#'
#' @description This function extracts the differential analysis results from
#'              the DETag object
#'
#' @param x An object of the DETag class
#' @param original A boolean indicating whether to return the original results.
#'                 If FALSE, the results will be converted to a data frame with
#'                 the following columns: logFC, pvalue, padj
#'                 (default: FALSE)
#'
#' @return A data frame containing the differential analysis results, with each
#'         row representing a gene and each column representing a feature
#' \itemize{
#' \item \code{row.names}: A character vector containing the gene names
#' \item \code{logFC}: A numeric vector containing the log2 fold change values
#' \item \code{pvalue}: A numeric vector containing the p-values
#' \item \code{padj}: A numeric vector containing the adjusted p-values
#' }
#'
#' @export
#'
#' @references
#' \insertRef{robinson2010edger}{IntegraTRN}
#'
#' \insertRef{love2014moderated}{IntegraTRN}
#'
#' @examples
#' # Use the package-provided example data
#' data("expMOList")
#' deTag <- expMOList$DERNAseq
#'
#' # Example 1: export the package default results
#' exportDE(deTag)
#'
#' # Example 2: export the original results as formatted by DESeq2 or edgeR
#' exportDE(deTag, original = TRUE)
#'
setGeneric("exportDE", function(x, original = FALSE) {
  standardGeneric("exportDE")
})
setMethod("exportDE", "DETag", function(x, original = FALSE) {
  # Validate the input
  if (!is.logical(original)) {
    stop("The original argument must be a boolean")
  } else {
    # Do nothing
  }
  # Extract the results
  if (original || x@method == ATAC_GRANGE) {
    # Differential accessibility analysis only results in a single format
    # so the original argument is ignored
    return(x@DEResult)
  } else {
    # Based on the type of the method, extract the results of count-based DE
    deResult <- x@DEResult
    if (x@method == DESEQ2) {
      deResult <- deResult %>% dplyr::select(DESEQ2_FIELDS)
    } else if (x@method == EDGER) {
      deResult <- deResult %>% dplyr::select(EDGER_FIELDS)
    } else {
      stop("The method is not supported") # never happen for a valid DETag obj
    }
    colnames(deResult) <- COUNT_DEFIELDS
    return(deResult)
  }
})


#' Export the normalized counts from a DETag object
#'
#' @aliases exportNormalizedCounts,DETag-method
#'
#' @description This function extracts the normalized counts from the DETag
#'              object, but this only applies to the count-based differential
#'              analysis
#'
#' @param x An object of the DETag class
#'
#' @return A matrix containing the normalized counts for the genes, with each
#'         row representing a gene and each column representing a sample
#'
#' @export
#'
#' @examples
#' # Use the package-provided example data
#' data("expMOList")
#' deTag <- expMOList$DERNAseq
#'
#' # Export the normalized counts
#' exportNormalizedCounts(deTag)
#'
setGeneric("exportNormalizedCounts", function(x) {
  standardGeneric("exportNormalizedCounts")
})
setMethod("exportNormalizedCounts", "DETag", function(x) {
  # Validate the input
  if (x@method != DESEQ2 && x@method != EDGER) {
    stop("The method is not supported")
  } else {
    # Do nothing
  }
  # Extract the normalized counts
  return(x@normalizedCounts)
})


# [END]
