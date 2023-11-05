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
#'
#' @inherit DETag
#'
#' @exportClass TOPTag
#'
#' @importFrom methods setClass
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
#' Reimand J, Isserlin R, Voisin V, et al. Pathway enrichment analysis and
#' visualization of omics data using g:Profiler, GSEA, Cytoscape and
#' EnrichmentMap. Nat Protoc. 2019;14(2):482-517. doi:10.1038/s41596-018-0103-9.
#'
#'
methods::setClass("TOPTag",
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
#'                10% of differential genes.
#'
#' @return A TOPTag object.
#'
#' @importFrom methods new
#' @importFrom dplyr filter arrange desc mutate %>%
#'
#' @export
#'
#' @examples
#' # Assuming that the object "deTag" is a DETag object
#'
#' # Example 1: Select the top 20% of differential genes with default cutoffs
#' \dontrun{
#' topTag <- TOPTag(deTag, topGenes = 0.2)
#' }
#'
#' # Example 2: Select the top 100 differential genes with default cutoffs
#' \dontrun{
#' topTag <- TOPTag(deTag, topGenes = 100)
#' }
#'
#' # Example 3: Select the top 20% of differential genes with custom cutoffs
#' \dontrun{
#' topTag <- TOPTag(deTag, logFCCutoff = 2, pCutoff = 0.01, topGenes = 0.2)
#' }
#'
TOPTag <- function(object, logFCCutoff = 1, pCutoff = 0.05, topGenes = 0.1) {
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

#' @rdname TOPTag-class
#'
#' @export
#'
#' @importFrom methods setMethod
#' @importFrom dplyr filter arrange desc %>%
#'
#' @method print TOPTag
#'
#' @description The print method for the TOPTag class.
#'
#' @inheritParams DETag-class
#'
#' @return NULL
#'
#' @examples
#' # Assuming that the object "topTag" is a TOPTag object
#' \dontrun{
#' print(topTag)
#' }
#'
#' # Or simply type the object name
#' \dontrun{
#' topTag
#' }
#'
methods::setMethod(
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
      head(5))
    cat("\nA snapshot of the top down-regulated genes:\n")
    print(object@DEResult %>%
      dplyr::filter(rank < 0) %>%
      dplyr::arrange(dplyr::desc(rank)) %>%
      head(5))
    cat("\n")
    cat("Method for differential expression analysis: ", object@method, "\n")
    cat("Log fold change cutoff: ", object@logFCCutoff, "\n")
    cat("Adjusted p-value cutoff: ", object@pCutoff)
  }
)


# [END]