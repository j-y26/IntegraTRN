# Purpose: Performing differential analysis on the omics data
# Author: Jielin Yang
# Date: 2023-10-30
# Version: 1.0
# Bugs and Issues: Currently, only DESeq2 is supported for differential analysis

#' DETag S4 class
#'
#' @keywords internal
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
#' @slot method A character string specifying the program used for the analysis
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

#' Constructor for the DETag class
#'
#' @keywords internal
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
  # Check the class of the input DEResult
  if (!inherits(DEResult, "data.frame")) {
    stop("The input DEResult must be a data frame")
  } else {
    # Do nothing
  }
  # Check the class of the input method
  if (!is.character(method)) {
    stop("The input method must be a character string")
  } else {
    # Do nothing
  }
  # Create the object
  newDETag <- new("DETag", DEResult = DEResult, method = method)

  # Return the object
  return(newDETag)
}

#### Need to define an S4 method to extract the differential results to
#### unified interface of dataframe
