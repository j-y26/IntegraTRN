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
edgeR_fields <- c("logFC", "logCPM", "PValue", "FDR", "Gene", "Annotation")
DESeq2_fields <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "GeneID", "Annotation")


#' @name DETag-class
#' @title DETag S4 class
#' @aliases DETag-class DETag
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
#'                results
#' @slot method A character string specifying the program used for the analysis,
#'              must be one of the valid METHODS defined
#'
#' @importFrom methods setClass
#'
#' @exportClass DETag
#'
#' @references
#' Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for
#' differential expression analysis of digital gene expression data.
#'
#' Love, M.I., Huber, W., and Anders, S. (2014). Moderated estimation of fold
#' change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 1–21.
#'
#' Advanced R by H. Wickham. Access: https://adv-r.hadley.nz/index.html
#'
methods::setClass("DETag", slots = c(
  DEResult = "data.frame",
  method = "character"
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
#' @examples
#' # Create an example data frame
#' deResult <- data.frame(
#'   gene = paste0("gene", seq_len(10)),
#'   logFC = runif(10),
#'   adj.P.Val = runif(10)
#' )
#'
#' # Create an object of the DETag class
#' deTag <- DETag(deResult, "DESeq2")
#'
#' # Validate the object
#' validateDETagSlots(deTag)
#'
validateDETagSlots <- function(DEResult, method) {
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
  return(invisible(NULL))
}


#' Constructor for the DETag class
#'
#' @description This function creates an object of the DETag class
#'
#' @param DEResult A data frame containing the differential expression analysis
#'                 results
#' @param method A character string specifying the program used for the analysis
#'
#' @return An object of the DETag class
#' \itemize{
#' \item \code{DEResult}: A data frame containing the differential expression
#'                        analysis results
#' \item \code{method}: A character string specifying the program used for the
#'                      analysis
#' }
#'
#' @importFrom methods new
#'
#' @export
#'
#' @references
#' Advanced R by H. Wickham. Access: https://adv-r.hadley.nz/index.html
#'
#' @examples
#' # Create an example data frame
#' deResult <- data.frame(
#'   gene = paste0("gene", seq_len(10)),
#'   logFC = runif(10),
#'   adj.P.Val = runif(10)
#' )
#'
#' # Create an object of the DETag class
#' deTag <- DETag(deResult, "DESeq2")
#'
#' # Check the class of the object
#' class(deTag)
#'
DETag <- function(DEResult, method) {
  # Validate the input
  validateDETagSlots(DEResult, method)
  # Create the object
  newDETag <- methods::new("DETag", DEResult = DEResult, method = method)

  # Return the object
  return(newDETag)
}


# Define the S4 method for the DETag class

#' @rdname DETag-class
#'
#' @method print DETag
#'
#' @description This function prints the DETag object
#'
#' @param x An object of the DETag class
#' @param ... Other arguments passed to the function
#'
#' @return NULL
#'
#' @importFrom methods standardGeneric setGeneric setMethod
#'
#' @export
#'
#' @references
#' Advanced R by H. Wickham. Access: https://adv-r.hadley.nz/index.html
#'
#' @examples
#' # Create an example data frame
#' deResult <- data.frame(
#'   gene = paste0("gene", seq_len(10)),
#'   logFC = runif(10),
#'   adj.P.Val = runif(10)
#' )
#'
#' # Create an object of the DETag class
#' deTag <- DETag(deResult, "DESeq2")
#'
#' # Print the object
#' print(deTag)
#'
print.DETag <- function(x, ...) {
  # Print the object
  cat("DETag object\n")
  cat("Method: ", x@method, "\n")
  cat("Differential analysis results:\n")
  print(exportDE(x))
  return(invisible(NULL))
}


# Generic function for extracting the differential analysis results
methods::setGeneric("exportDE", function(x, original = FALSE) {
  standardGeneric("exportDE")
})

#' @rdname DETag-class
#'
#' @method exportDE DETag
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
#' @importFrom methods standardGeneric setGeneric setMethod
#'
#' @export
#'
#' @references
#' Advanced R by H. Wickham. Access: https://adv-r.hadley.nz/index.html
#'
#' Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for
#' differential expression analysis of digital gene expression data.
#' Bioinformatics. 2010 Jan 1;26(1):139-40. doi: 10.1093/bioinformatics/btp616.
#' Epub 2009 Nov 11. PMID: 19910308; PMCID: PMC2796818.
#'
#' Love, M.I., Huber, W., and Anders, S. (2014). Moderated estimation of fold
#' change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 1–21.
#'
#' @examples
#' # Example 1: export the package default results
#'
#' \dontrun{
#' # Assuming the deTag object is already created from a differential analysis
#'
#' # Export the results
#' exportDE(deTag)
#' }
#'
#' # Example 2: export the original results
#'
#' \dontrun{
#' # Export the results while keeping the original results from DESeq2 or edgeR
#' exportDE(deTag, original = TRUE)
#' }
#'
methods::setMethod("exportDE", "DETag", function(x, original = FALSE) {
  # Validate the input
  if (!is.logical(original)) {
    stop("The original argument must be a boolean")
  } else {
    # Do nothing
  }
  # Extract the results
  if (original) {
    return(x@DEResult)
  } else {
    # Based on the type of the method, extract the results
    if (x@method %in% c(DESEQ2, EDGER)) {
      # DESeq2 or edgeR
      # Extract the results
      if (x@method == DESEQ2) {
        # DESeq2
        # Extract the results
        deResult <- x@DEResult[, DESeq2_fields]
        # Rename the columns
        colnames(deResult) <- c("logFC", "pvalue", "padj")
      } else {
        # edgeR
        # Extract the results
        deResult <- x@DEResult[, edgeR_fields]
        # Rename the columns
        colnames(deResult) <- c("logFC", "logCPM", "pvalue", "padj")
      }
      # Return the results
      return(deResult)
    } else {
      # ATAC-GRanges
      # Extract the results
      deResult <- x@DEResult
      # Rename the columns
      colnames(deResult) <- c("logFC", "pvalue", "padj")
      # Return the results
      return(deResult)
    }
  }
})
