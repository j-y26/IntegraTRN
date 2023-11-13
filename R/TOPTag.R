# Purpose: TOPTag S4 class for storing top differential genes
# Author: Jielin Yang
# Date: 2023-11-04
# Version: 1.0
# Bugs and Issues: None


#' @title TOPTag S4 class for storing top differential genes
#'
#' @description The TOPTag S4 class inherits the DETag class and is a more
#'              specific class for storing top differential genes. It contains
#'              all slots of the DETag class, with stringent filtering on the
#'              genes involved. Additional slots are added as a more convenient
#'              way for users to access the top differential genes.
#'
#'
#' @slot DEResult A data frame containing the differential expression results
#'                for only the top differential genes. The columns resemble the
#'                predefined format in the DETag class. An additional column
#'                "rank" is added to indicate the rank of the gene in the
#'                differential expression results.
#' \itemize{
#'  \item {logFC}{The log fold change of the gene}
#'  \item {pvalue}{The p-value of the gene}
#'  \item {padj}{The adjusted p-value of the gene}
#'  \item {piValue}{The pi-value of the gene, calculated by the formula
#'                 -log10(padj) * sign(logFC)}
#'  \item {rank}{The rank of the gene in the differential expression results.
#'               Positive for up-regulated genes and negative for down-
#'               regulated genes. A rank of 1 indicates the most up-regulated
#'               gene and a rank of -1 indicates the most down-regulated gene.}
#' }
#' @slot method A character string indicating the method used to calculate the
#'              differential expression results.
#' @slot normalizedCounts A matrix containing the normalized counts for
#'                        only the top differential genes, consistent with the
#'                        the DEResult slot.
#' @slot logFCCutoff A numeric value indicating the absolute log2 fold change
#'                   cutoff used to filter the differential expression results.
#' @slot pCutoff A numeric value indicating the adjusted p-value cutoff used to
#'               filter the differential expression results.
#' @slot topGenes A numeric value indicating either the fraction (0-1) of top
#'                differential genes or the number (1-Inf) of top differential
#'                genes. If the number specified is greater than the number of
#'                DE genes based on logFCCutoff and pCutoff, then topGenes will
#'                be set to the number of DE genes.
#'
#' @exportClass TOPTag
#'
#' @references
#' \insertRef{love2014moderated}{IntegraTRN}
#'
#' \insertRef{robinson2010edger}{IntegraTRN}
#'
#' \insertRef{reimand2019pathway}{IntegraTRN}
#'
setClass("TOPTag",
  contains = "DETag",
  slots = c(
    DEResult = "data.frame",
    method = "character",
    normalizedCounts = "matrix",
    logFCCutoff = "numeric",
    pCutoff = "numeric",
    topGenes = "numeric"
  )
)


#' Constructor for TOPTag class
#'
#' @aliases TOPTag
#'
#' @description The constructor for the TOPTag class.
#'
#' @param object A DETag object.
#' @param logFCCutoff A numeric value indicating the absolute log2 fold change
#'                    cutoff used to filter the differential expression results.
#' @param pCutoff A numeric value indicating the adjusted p-value cutoff used to
#'                filter the differential expression results.
#' @param topGenes A numeric value indicating either the fraction (0-1) of top
#'                differential genes or the number (1-Inf) of top differential
#'                genes. If the number specified is greater than the number of
#'                DE genes based on logFCCutoff and pCutoff, then topGenes will
#'                be set to the number of DE genes. Default to select the top
#'                10\% of differential genes.
#' @param direction A character string indicating the direction of the
#'                  differential expression. Default to "both". Other options
#'                  include "up" and "down".
#'
#' @return A TOPTag object.
#'
#' @importFrom dplyr filter arrange desc mutate %>%
#'
#' @export
#'
#' @references
#' \insertRef{dplyr}{IntegraTRN}
#'
#' \insertRef{reimand2019pathway}{IntegraTRN}
#'
#' @examples
#' # Use the package provided example data
#' data("expMOList")
#' # Extract an DETag object from RNA-seq data
#' deTag <- expMOList$DERNAseq
#'
#' # Example 1: Select the top 20% of differential genes with default cutoffs
#' topTag <- TOPTag(deTag, topGenes = 0.2)
#' topTag
#'
#' # Example 2: Select the top 50 differential genes with default cutoffs
#' \dontrun{
#' topTag <- TOPTag(deTag, topGenes = 50)
#' topTag
#' }
#'
#' # Example 3: Select the top 20% of differential genes with custom cutoffs
#' topTag <- TOPTag(deTag, logFCCutoff = 0, pCutoff = 0.01, topGenes = 0.2)
#' topTag
#'
#' # Example 4: Select the top 20% of up-regulated genes with default cutoffs
#' topTag <- TOPTag(deTag, topGenes = 0.2, direction = "up")
#' topTag
#'
TOPTag <- function(object,
                   logFCCutoff = 1,
                   pCutoff = 0.05,
                   topGenes = 0.1,
                   direction = "both") {
  # Validating the inputs
  if (!is(object, "DETag")) {
    stop("The input object must be a DETag object.")
  } else if (!all(is.numeric(c(logFCCutoff, pCutoff, topGenes)))) {
    stop("The logFCCutoff, pCutoff and topGenes must be numeric.")
  } else if (topGenes < 0) {
    stop("The topGenes must be a positive number.")
  } else {
    # Do nothing
  }

  # logFCCutoff will be coerced to a positive number
  logFCCutoff <- abs(logFCCutoff)
  # If topGenes > 1, it must be an integer
  if (topGenes > 1) {
    topGenes <- round(topGenes)
  } else {
    # Do nothing
  }

  # Filter by logFCCutoff and pCutoff only
  deGenes <- object %>%
    exportDE() %>%
    dplyr::filter(
      abs(logFC) >= logFCCutoff,
      padj <= pCutoff
    )
  # Directional filtering
  if (direction == "up") {
    deGenes <- deGenes %>% dplyr::filter(logFC > 0)
  } else if (direction == "down") {
    deGenes <- deGenes %>% dplyr::filter(logFC < 0)
  } else {
    # Do nothing
  }

  if (topGenes > 1 && topGenes > nrow(deGenes)) {
    warning("Selected number of top genes is greater than the number of DE
             genes. Selecting all DE genes.")
    topGenes <- nrow(deGenes)
  } else {
    # Do nothing
  }

  # Rank the genes
  deGenes <- deGenes %>%
    dplyr::mutate(piValue = -log(padj, 10) * sign(logFC)) %>%
    dplyr::arrange(dplyr::desc(piValue)) %>%
    dplyr::mutate(rank = ifelse(logFC > 0,
      seq(1, nrow(.[logFC > 0, ])),
      seq(-nrow(.[logFC < 0, ]), -1)
    ))
  # Select the top genes
  if (topGenes > 1) {
    # Proportional selection for up and down regulated genes
    numTopUp <- round(topGenes * (sum(deGenes$logFC > 0) / nrow(deGenes)))
    numTopDown <- topGenes - numTopUp
    topDE <- rbind(
      deGenes %>%
        dplyr::filter(rank > 0) %>%
        dplyr::filter(rank <= numTopUp),
      deGenes %>%
        dplyr::filter(rank < 0) %>%
        dplyr::filter(rank >= -numTopDown)
    )
  } else {
    # Select the top genes based on the fraction
    topDE <- deGenes %>%
      dplyr::filter(abs(rank) <= round(nrow(deGenes) * topGenes))
  }

  # Construct the TOPTag object
  topTag <- new("TOPTag",
    DEResult = topDE,
    method = object@method,
    normalizedCounts = object@normalizedCounts[rownames(topDE), ],
    logFCCutoff = logFCCutoff,
    pCutoff = pCutoff,
    topGenes = topGenes
  )

  # Clean up
  rm(deGenes, topDE)

  return(topTag)
}


# Define the S4 method for the TOPTag class

#' print TOPTag
#'
#' @aliases print,TOPTag-method
#'
#' @export
#'
#' @importFrom dplyr filter arrange desc %>%
#'
#' @method print TOPTag
#'
#' @description The print method for the TOPTag class.
#'
#' @param object A TOPTag object.
#'
#' @return NULL
#'
#' @examples
#' # Use the package provided example data
#' data("expMOList")
#'
#' # Extract an DETag object from RNA-seq data
#' deTag <- expMOList$DERNAseq
#'
#' # Create a TOPTag object
#' topTag <- TOPTag(deTag, topGenes = 0.2)
#'
#' # Print the TOPTag object
#' print(topTag)
#'
#' # or simply type the object name
#' topTag
#'
setMethod(
  "show", "TOPTag",
  function(object) {
    if (object@topGenes > 1) {
      cat("A TOPTag object.\n")
      cat("The top", object@topGenes, "differential genes are selected.\n\n")
    } else {
      cat(
        "A TOPTag object with", nrow(object@DEResult),
        "top differential genes.\n"
      )
      cat(
        "The top", object@topGenes * 100,
        "% differential genes are selected.\n\n"
      )
    }
    cat("A snapshot of the top up-regulated genes:\n")
    print(object@DEResult %>%
      dplyr::filter(rank > 0) %>%
      dplyr::arrange(rank) %>%
      utils::head(5))
    cat("\nA snapshot of the top down-regulated genes:\n")
    print(object@DEResult %>%
      dplyr::filter(rank < 0) %>%
      dplyr::arrange(dplyr::desc(rank)) %>%
      utils::head(5))
    cat("\n")
    cat("Method for differential expression analysis: ", object@method, "\n")
    cat("Log fold change cutoff: ", object@logFCCutoff, "\n")
    cat("Adjusted p-value cutoff: ", object@pCutoff)
  }
)


#' @rdname TOPTag-class
#'
#' @keywords internal
#'
#' @title Filter genes by a names vector
#'
#' @description This function filters the top differential genes by a names
#'              vector
#'
#' @param object A TOPTag object.
#' @param names A character vector containing the names of the genes to be
#'              filtered.
#'
#' @return A TOPTag object
#'
#' @examples
#' # Use the package provided example data
#' data("expMOList")
#'
#' # Extract an DETag object from RNA-seq data
#' deTag <- expMOList$DERNAseq
#'
#' # Create a TOPTag object
#' topTag <- TOPTag(deTag, logFCCutoff = 0, topGenes = 0.2)
#'
#' # Filter the top differential genes by a names vector
#' topTag <- filterGenes(topTag, c("B2M", "VSIR", "HAND2-AS1"))
#'
setGeneric("filterGenes", function(object, names) {
  standardGeneric("filterGenes")
})
setMethod("filterGenes", "TOPTag", function(object, names) {
  object@DEResult <- object@DEResult[rownames(object@DEResult) %in% names, ]
  object@normalizedCounts <- object@normalizedCounts[
    rownames(object@DEResult),
  ]
  return(object)
})


#' Export the differential analysis results from a TopTag object
#'
#' @aliases exportDE,TOPTag-method
#'
#' @description This function exports the differential analysis results from a
#'              TopTag object. Overwrites the exportDE function in the DETag
#'              class due to different nomenclature of the columns.
#'
#' @inheritParams exportDE
#'
#' @return A data frame containing the differential analysis results.
#'
#' @export
#'
#' @examples
#' # Use the package provided example data
#' data("expMOList")
#'
#' # Extract an DETag object from RNA-seq data
#' deTag <- expMOList$DERNAseq
#'
#' # Create a TOPTag object
#' topTag <- TOPTag(deTag, logFCCutoff = 0, topGenes = 0.2)
#'
#' # Export the differential analysis results
#' exportDE(topTag)
#'
setMethod("exportDE", "TOPTag", function(x, original = FALSE) {
  # Ignore the original argument
  return(x@DEResult)
})


# [END]
