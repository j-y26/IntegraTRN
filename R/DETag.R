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
  newDETag <- new("DETag", DEResult = DEResult, method = method)

  # Return the object
  return(newDETag)
}

#### Need to define an S4 method to extract the differential results to
#### unified interface of dataframe

exportDE