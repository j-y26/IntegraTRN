# Purpose: Constructing transcriptional regulatory network
# Author: Jielin Yang
# Date: 2023-11-07
# Version: 1.0
# Bugs and Issues: None


# Define global variables
NETWORK_FIELD <- c("regulator", "target", "regulatorType")


#' Validate external interaction data
#'
#' @keywords internal
#'
#' @description This function validates the external interaction data provided
#'              by the user and returns a list of validated data.
#'
#' @param adjList A list of interactions in adjacent list format
#'
#' @return A list of validated interactions in adjacent list format
#'
validateInteractionAdjList <- function(adjList) {
  if (!is.list(adjList) || length(adjList) < 2) {
    stop("The input must be a list of at least two elements.")
  } else {
    # Continue
  }

  if (!all(NETWORK_FIELD[1:2] %in% names(adjList))) {
    warning(paste0(
      "Invalid element names provided. Setting the first element ",
      "as regulator and second as target."
    ))
    names(adjList)[1:2] <- NETWORK_FIELD[1:2]
  } else {
    # Do nothing
  }
  return(adjList[NETWORK_FIELD[1:2]])
}


#' Loading external interaction data into the MOList object
#'
#' @aliases loadExtInteractions
#'
#' @description This function loads external interaction data into the MOList
#'              object. The external data can be curated non-tissue/sample
#'              specific interaction data, or any interactions that the users
#'              has determined to encompass a global picture of the regulatory
#'              interactions of the target genes.
#'
#' @param objMOList A MOList object containing the omics data
#' @param miR2Genes A list containing the target genes and their regulatory
#'                  miRNAs. Must contain elements "regulator" for the miRNAs
#'                  and "target" for the target genes. The two elements must
#'                  be character vectors of the same length.
#' \itemize{
#'  \item \code{regulator}: A character vector containing the miRNAs
#'  \item \code{target}: A character vector containing the target genes
#' }
#' @param tf2Genes A list containing the target genes and their regulatory
#'                 transcription factors. Must contain elements "regulator" for
#'                 the transcription factors and "target" for the target genes.
#'                 The two elements must be character vectors of the same
#'                 length. See the examples for detailed format.
#'
#' @return An object of class MOList, with a element "extInteractions" added,
#'         which is a list of lists. Each lower level list must follow the
#'         format of the input data, i.e., must contain elements "regulator"
#'         and "target" in the adjacent list format. The names of the regulator
#'         and target elements should correspond to the names provided with the
#'         input data. See the examples for details.
#' \itemize{
#' \item \code{miR2Genes}: A list containing the target genes
#'                              and their regulatory miRNAs
#' \item \code{tf2Genes}: A list containing the target genes
#'                             and their regulatory transcription factors
#' }
#'
#' @export
#'
#' @examples
#' # Create a sample RNAseq data
#' RNAseq <- matrix(sample(1:1000, 1000, replace = TRUE), ncol = 10)
#' colnames(RNAseq) <- paste0("sample_", seq_len(ncol(RNAseq)))
#' rownames(RNAseq) <- paste0("gene_", seq_len(nrow(RNAseq)))
#' RNAGroupBy <- rep(c("A", "B"), each = 5)
#'
#' # Create a myMOList object
#' myMOList <- MOList(RNAseq = RNAseq, RNAGroupBy = RNAGroupBy)
#'
#' # Generate some example interaction data
#' miR2Genes <- list(
#'   regulator = c("miR_1", "miR_2", "miR_3", "miR_4", "miR_5"),
#'   target = c("gene_1", "gene_2", "gene_3", "gene_4", "gene_5")
#' )
#' tf2Genes <- list(
#'   regulator = c("TF_1", "TF_2", "TF_3", "TF_4", "TF_5"),
#'   target = c("gene_1", "gene_2", "gene_3", "gene_4", "gene_5")
#' )
#' # Load the external interaction data into the MOList object
#' myMOList <- loadExtInteractions(
#'   myMOList,
#'   miR2Genes = miR2Genes,
#'   tf2Genes = tf2Genes
#' )
#'
loadExtInteractions <- function(objMOList,
                                miR2Genes = NULL,
                                tf2Genes = NULL) {
  # Check the input data
  if (is.null(miR2Genes) && is.null(tf2Genes)) {
    stop("Please provide at least one type of the interaction data.")
  } else {
    # Continue
  }

  if (!is.null(miR2Genes)) {
    # In case where the names of the elements of the list are not provided
    # correctly, set the first element to be the regulator and the second
    # element to be the target
    miR2Genes <- validateInteractionAdjList(miR2Genes)
  } else {
    # Do nothing
  }
  if (!is.null(tf2Genes)) {
    tf2Genes <- validateInteractionAdjList(tf2Genes)
  } else {
    # Do nothing
  }

  # Setting the external interactions to the MOList object
  objMOList$extInteractions <- list(
    miR2Genes = miR2Genes,
    tf2Genes = tf2Genes
  )

  return(objMOList)
}


#' Setting cutoffs for the omics data
#'
#' @aliases setOmicCutoffs
#'
#' @description This function sets the cutoffs for the omics data. The cutoffs
#'              are used to determine the key differentially regulated/expressed
#'              genes/miRNAs/TFs. The users can set the cutoffs based on their
#'              interpretation during the exploratory analysis on the
#'              differential analysis process. Here, the cutoffs are set to
#'              define a set of genes/miRNAs/TFs that are the key to such
#'              differential regulation/expression, and only these genes/miRNAs/
#'              TFs will be used to construct the network. Any values set on
#'              data that is not available will be ignored.
#'
#' @param rnaAdjPval The adjusted p-value cutoff for the RNAseq data, must be
#'                   a numeric value between 0 and 1. Default is 0.05.
#' @param rnaLogFC The log fold change cutoff for the RNAseq data, must be a
#'                 positive number. Default is 0.
#' @param rnaTopGenes A numeric value indicating either the fraction (0-1) of
#'                    top differential genes or the number (1-Inf) of top
#'                    differential genes. If the number specified is greater
#'                    than the number of DE genes based on logFCCutoff and
#'                    pCutoff, then topGenes will be set to the number of DE
#'                    genes. This is to select the top differential mRNAs.
#'                    To include all the differential mRNAs, set this to 1.
#'                    Default is 1.
#' @param smallRNAAdjPval The adjusted p-value cutoff for the smallRNAseq data.
#'                        Must be a numeric value between 0 and 1. Default is
#'                        0.05.
#' @param smallRNALogFC The log fold change cutoff for the smallRNAseq data.
#'                      Must be a positive number. Default is 0.
#' @param smallRNATopGenes A numeric value indicating either the fraction (0-1)
#'                         of top differential genes or the number (1-Inf) of
#'                         top differential genes. If the number specified is
#'                         greater than the number of DE genes based on
#'                         logFCCutoff and pCutoff, then topGenes will be set
#'                         to the number of DE genes. This is to select the top
#'                         differential small RNAs. To include all the
#'                         differential small RNAs, set this to 1. Default is
#'                         1.
#' @param proteomicsAdjPval The adjusted p-value cutoff for the proteomics data.
#'                          Must be a numeric value between 0 and 1. Default is
#'                          0.05.
#' @param proteomicsLogFC The log fold change cutoff for the proteomics data.
#'                        Must be a positive number. Default is 0.
#' @param atacMotifAdjPval The adjusted p-value cutoff for the ATACseq motif
#'                         enrichment analysis. Must be a numeric value between
#'                         0 and 1. Default is 0.05.
#' @param atacMotifPval The p-value cutoff for the ATACseq motif enrichment
#'                      analysis. This value is defaulted to NULL, but if the
#'                      user specifies this value, then the adjusted p-value
#'                      cutoff will be ignored.
#' @param atacMotifLogFC The log2 fold enrichment cutoff for the ATACseq motif
#'                       enrichment analysis.
#'
#' @return An OMICutoffs object containing the cutoffs for the omics data. This
#'         is essentially a list with defined elements. Once setOmicCutoffs is
#'         called, users can freely access or modify the cutoffs by using the
#'         $ operator. See the examples for details.
#'
#' @export
#'
#' @examples
#' # Example 1: Set the cutoffs for the omics data
#' omiCutoffs <- setOmicCutoffs(
#'   rnaAdjPval = 0.05,
#'   rnaLogFC = 2,
#'   rnaTopGenes = 0.1, # Top 10% of DE genes
#'   smallRNAAdjPval = 0.05,
#'   smallRNALogFC = 2,
#'   smallRNATopGenes = 200, # Top 200 DE small RNAs
#'   proteomicsAdjPval = 0.05,
#'   proteomicsLogFC = 1
#' )
#'
#' # Example 2: Access the cutoffs for the RNAseq data
#' omiCutoffs$rnaAdjPval
#' omiCutoffs$rnaLogFC
#'
#' # Example 3: Modify the cutoffs for the proteomics data
#' omiCutoffs$proteomicsAdjPval <- 0.01
#' omiCutoffs$proteomicsLogFC <- 1.5
#'
setOmicCutoffs <- function(rnaAdjPval = 0.05,
                           rnaLogFC = 0,
                           rnaTopGenes = 1,
                           smallRNAAdjPval = 0.05,
                           smallRNALogFC = 0,
                           smallRNATopGenes = 1,
                           proteomicsAdjPval = 0.05,
                           proteomicsLogFC = 0,
                           atacMotifAdjPval = 0.05,
                           atacMotifPval = NULL,
                           atacMotifLogFC = NULL) {
  # Check the inputs
  inputs <- c(
    rnaAdjPval, rnaLogFC, rnaTopGenes,
    smallRNAAdjPval, smallRNALogFC, smallRNATopGenes,
    proteomicsAdjPval, proteomicsLogFC,
    atacMotifAdjPval
  )
  if (!all(is.null(inputs) | is.numeric(inputs))) {
    stop("The cutoffs must be numeric.")
  } else if (!all(is.null(inputs) | inputs >= 0)) {
    stop("The cutoffs must be a positive number.")
  } else {
    # Continue
  }

  # Setting the cutoffs
  omiCutoffs <- list(
    rnaAdjPval = rnaAdjPval,
    rnaLogFC = rnaLogFC,
    rnaTopGenes = rnaTopGenes,
    smallRNAAdjPval = smallRNAAdjPval,
    smallRNALogFC = smallRNALogFC,
    smallRNATopGenes = smallRNATopGenes,
    proteomicsAdjPval = proteomicsAdjPval,
    proteomicsLogFC = proteomicsLogFC,
    atacMotifAdjPval = atacMotifAdjPval,
    atacMotifPval = atacMotifPval,
    atacMotifLogFC = atacMotifLogFC
  )

  # Setting the class, an internal S3 class
  class(omiCutoffs) <- "OMICutoffs"

  return(omiCutoffs)
}


#' Print method for the OMICutoffs object
#'
#' @aliases print.OMICutoffs
#'
#' @description This function prints the cutoffs for the omics data to the
#'              console.
#'
#' @param x An OMICutoffs object
#' @param ... Other arguments passed to the print function
#'
#' @return The cutoffs for the omics data
#'
#' @export
#'
#' @examples
#' omiCutoffs <- setOmicCutoffs()
#' print(omiCutoffs)
#'
print.OMICutoffs <- function(x, ...) {
  cat("RNAseq adjusted p-value cutoff: <", x$rnaAdjPval, "\n")
  cat("RNAseq log fold change cutoff: >", x$rnaLogFC, "\n")
  if (x$rnaTopGenes <= 1) {
    cat("RNAseq selecting top", x$rnaTopGenes * 100, "% of DE genes", "\n")
  } else {
    cat("RNAseq selecting top", x$rnaTopGenes, "DE genes", "\n")
  }
  if (!is.null(x$smallRNAAdjPval)) {
    cat("smallRNAseq adjusted p-value cutoff: <", x$smallRNAAdjPval, "\n")
    cat("smallRNAseq log fold change cutoff: >", x$smallRNALogFC, "\n")
    if (x$smallRNATopGenes <= 1) {
      cat(
        "smallRNAseq selecting top", x$smallRNATopGenes * 100, "% of DE genes",
        "\n"
      )
    } else {
      cat("smallRNAseq selecting top", x$smallRNATopGenes, "DE genes", "\n")
    }
  } else {
    cat("No smallRNAseq value set.", "\n")
  }
  if (!is.null(x$proteomicsAdjPval)) {
    cat("Proteomics adjusted p-value cutoff: <", x$proteomicsAdjPval, "\n")
    cat("Proteomics log fold change cutoff: >", x$proteomicsLogFC, "\n")
  } else {
    cat("No proteomics value set.", "\n")
  }
  if (is.null(x$atacMotifPval)) {
    cat("ATACseq motif adjusted p-value cutoff: <", x$atacMotifAdjPval, "\n")
  } else {
    cat("ATACseq motif p-value cutoff: <", x$atacMotifPval, "\n")
  }
  atacFE <- ifelse(is.null(x$atacMotifLogFC), 0, x$atacMotifLogFC)
  cat("ATACseq motif log fold enrichment cutoff: >", atacFE, "\n")
}


#' Select miRNA inverse correlation
#'
#' @description This function selects the miRNAs in which its expression is
#'              inversely correlated with the expression of the target gene.
#'              Regulators that are not miRNAs will be retained.
#'
#' @param exprAdjList A list representing the regulator-target interactions
#' \itemize{
#' \item \code{regulator}: A character vector containing the regulators
#' \item \code{target}: A character vector containing the targets
#' }
#' @param deResultRNA A data frame containing the differential RNAseq data
#' @param deResultSmallRNA A data frame containing the differential smallRNAseq
#' @param sncAnno A list of annotations for small RNAs
#' @param smallRNATypes A character vector containing the small RNA types that
#'                      the user wants to use to construct the network. The
#'                      available types are "miRNA", "piRNA", "snRNA", "snoRNA",
#'                      "tRNA", and "circRNA".
#'
#' @return A list with filtered regulator-target interactions
#'
#' @importFrom dplyr filter %>%
#'
filtermiRNAinverseCorr <- function(exprAdjList,
                                   deResultRNA,
                                   deResultSmallRNA,
                                   sncAnno,
                                   smallRNATypes = SMALLRNA_CATEGORIES) {
  if ("miRNA" %in% smallRNATypes) {
    # Extract only inverse small RNA - mRNA interactions
    # i.e., downregulated small RNA - upregulated mRNA
    genesmiRNA <- extractDirectionalGenes(deResultSmallRNA) %>%
      lapply(intersect, sncAnno[["miRNA"]])
    genesRNA <- extractDirectionalGenes(deResultRNA)

    downregmiRNA <- genesmiRNA$down
    upregmiRNA <- genesmiRNA$up
    downregRNA <- genesRNA$down
    upregRNA <- genesRNA$up

    rm(genesmiRNA, genesRNA)

    # Extract the predicted interactions
    inverseRelation <- (exprAdjList$regulator %in% downregmiRNA &
      exprAdjList$target %in% upregRNA) |
      (exprAdjList$regulator %in% upregmiRNA &
        exprAdjList$target %in% downregRNA)
    notmiRNA <- !(exprAdjList$regulator %in% sncAnno[["miRNA"]])
    exprAdjList <- lapply(
      exprAdjList,
      function(x) x[inverseRelation | notmiRNA]
    )
  } else {
    # Do nothing
  }
  return(exprAdjList)
}


#' Predict small RNA - mRNA interactions
#'
#' @keywords internal
#'
#' @description This function predicts the small RNA - mRNA interactions based
#'              on the smallRNAseq data and the RNAseq data.
#'
#' @param objMOList A MOList object containing the omics data
#' @param omiCutoffs A OMICutoffs object containing the cutoffs for the omics
#'                   data
#' @param rnaTopTag A TOPTag object containing the top differential genes from
#'                  the RNAseq data
#' @param smallRNATypes A character vector containing the small RNA types that
#'                      the user wants to use to predict the interactions. The
#'                      available types are "miRNA", "piRNA", "snRNA", "snoRNA",
#'                      "tRNA", and "circRNA".
#' @param ntree The number of trees to grow in the random forest model
#' @param nthreads The number of threads to use for parallel processing
#' @param treeMethod The method to use for the random forest model. See the
#'                  documentation for the randomForest package for details.
#' @param seed The seed to use for the random forest model
#'
#' @return a list of predicted interactions combined for all specified small
#'         RNA types as the regulators
#'
predictSncmRNAInteractions <- function(objMOList, omiCutoffs, rnaTopTag,
                                       smallRNATypes = SMALLRNA_CATEGORIES,
                                       ntree = 1000, nthreads = 1,
                                       treeMethod = "RF", seed = 91711) {
  if (is.null(objMOList$matchingRNAsmallRNA)) {
    stop("No sample matching between the RNAseq and smallRNAseq data. Please
    perform sample matching before constructing the network.")
  } else {
    # Continue
  }

  # Generates a TOPTag object for the small RNAs
  smallRNATopTag <- TOPTag(objMOList$DEsmallRNAseq,
    logFCCutoff = omiCutoffs$smallRNALogFC,
    pCutoff = omiCutoffs$smallRNAAdjPval,
    topGenes = omiCutoffs$smallRNATopGenes,
    direction = "both"
  )

  # Define annotations for small RNA
  if (all(objMOList$annoSncRNA == HUMAN)) {
    sncAnno <- SNCANNOLIST_HSAPIENS
  } else {
    sncAnno <- objMOList$annoSncRNA
  }

  # Predict the interactions for each small RNA type
  predInteract <- predictSmallRNAmRNAcoExpr(
    mRNATopTag = rnaTopTag,
    smallRNATopTag = smallRNATopTag,
    smallRNATypes = smallRNATypes,
    annoSncRNA = sncAnno,
    matchingRNAsmallRNA = objMOList$matchingRNAsmallRNA,
    ntree = ntree,
    nthreads = nthreads,
    treeMethod = treeMethod,
    seed = seed
  )
  return(predInteract)
}


#' Intersect interaction lists
#'
#' @description This function intersects two lists of interactions in
#'              adjacent list format.
#'
#' @note At least one of the lists must be non-empty.
#'
#' @param adjList1 A list of interactions in adjacent list format
#' \itemize{
#' \item \code{regulator}: A character vector containing the regulators
#' \item \code{target}: A character vector containing the targets
#' }
#' @param adjList2 A list of interactions in adjacent list format
#'
#' @return A list of interactions in adjacent list format, with the interactions
#'         that are common to both lists
#'
intersectInteractions <- function(adjList1, adjList2) {
  # Simply returns the other list if one of the list is empty
  if (is.null(adjList1)) {
    return(adjList2)
  } else if (is.null(adjList2)) {
    return(adjList1)
  } else {
    # Continue
  }

  # Paired interactions only, so define a new element for easy comparison
  adjList1$pair <- paste0(adjList1$regulator, "_", adjList1$target)
  adjList2$pair <- paste0(adjList2$regulator, "_", adjList2$target)

  # Intersect the interactions
  commonInteract <- intersect(adjList1$pair, adjList2$pair)
  sel <- adjList1$pair %in% commonInteract

  intersectedList <- list(
    regulator = adjList1$regulator[sel],
    target = adjList1$target[sel]
  )

  return(intersectedList)
}


#' Filter target genes by proteomics data
#'
#' @keywords internal
#'
#' @description This function filters the target genes by proteomics data.
#'              Only genes that shows consistent expression changes in both
#'              RNAseq and proteomics data will be used to construct the
#'              network. This is to ensure high confidence in the target genes.
#'
#' @param rnaTopTag A TOPTag object containing the top differential genes from
#'                  the RNAseq data
#' @param protTopTag A TOPTag object containing the top differential genes from
#'                   the proteomics data
#' @param mapping A data frame containing the mapping between genes and proteins
#' \itemize{
#' \item \code{gene}: A character vector containing the gene names
#' \item \code{protein}: A character vector containing the protein names
#' }
#'
#' @return A TOPTag object containing the top differential genes from the
#'         RNAseq data, with the target genes filtered by proteomics data
#'
#' @importFrom dplyr %>%
#'
filterTargetGenes <- function(rnaTopTag, protTopTag, mapping) {
  # Filter the target genes by proteomics data
  validProteins <- protTopTag %>%
    exportDE() %>%
    rownames()
  validTargets <- mapping$gene[mapping$protein %in% validProteins]

  # Redefine the target genes
  rnaTopTag <- filterGenes(rnaTopTag, validTargets)
  return(rnaTopTag)
}


#' Combine small RNA - mRNA interactions from different sources
#'
#' @keywords internal
#'
#' @description This function combines the small RNA - mRNA interactions from
#'              different sources, including the predicted interactions and the
#'              curated interactions. The function also annotates the type of
#'              the regulators. For miRNA only, the functions also filters the
#'              interactions by inverse correlation.
#'
#' @param predInteract A list of predicted interactions in adjacent list format
#' \itemize{
#' \item \code{regulator}: A character vector containing the regulators
#' \item \code{target}: A character vector containing the targets
#' }
#'
#' @param extInteract A list of curated interactions in adjacent list format
#' @param sncAnno A list of annotations for small RNAs
#' @param deResultRNA A data frame containing the differential RNAseq data
#' @param deResultSmallRNA A data frame containing the differential smallRNAseq
#'
#' @return A list of interactions in adjacent list format, with the interactions
#'         that are common to both lists, with the regulator type annotated,
#'         and with the miRNA - mRNA interactions filtered. If only one type
#'         of interactions is provided, then the function will return the
#'         interactions directly.
#'
combineSncInteractions <- function(predInteract,
                                   extInteract,
                                   sncAnno,
                                   deResultRNA,
                                   deResultSmallRNA,
                                   smallRNATypes = SMALLRNA_CATEGORIES) {
  predicted <- is.null(predInteract)
  external <- is.null(extInteract)

  if (predicted && external) {
    return(NULL)
  } else if (predicted) {
    # No predicted interactions
    sncInteract <- extInteract
  } else if (external) {
    # No curated interactions
    sncInteract <- predInteract
  } else {
    # Identify miRNAs to intersect, while keeping the non-miRNA regulations
    selmiRNA <- predInteract$regulator %in% sncAnno[["miRNA"]]
    predmiRNA <- lapply(predInteract, function(x) x[selmiRNA])
    predNonmiRNA <- lapply(predInteract, function(x) x[!selmiRNA])
    # Intersect the miRNA - mRNA interactions
    miRNAmRNA <- intersectInteractions(predmiRNA, extInteract)
    # Then restore the list of interactions
    sncInteract <- mapply(c, predNonmiRNA, miRNAmRNA, SIMPLIFY = FALSE)
  }
  # Executed only if some interactions are available and has small RNA expr
  # Filter for inverse correlation on miRNA - mRNA interactions
  if (!is.null(deResultSmallRNA)) {
    sncInteract <- filtermiRNAinverseCorr(
      exprAdjList = sncInteract,
      deResultRNA = deResultRNA,
      sncAnno = sncAnno,
      deResultSmallRNA = deResultSmallRNA,
      smallRNATypes = smallRNATypes
    )
  } else {
    # Filter interactions simply by targets
  }
  # Annotate the small RNA types
  if (!is.null(deResultSmallRNA)) {
    sncInteract$regulatorType <- findGeneType(sncInteract$regulator, sncAnno)
  } else if (predicted && !external) {
    sncInteract$regulatorType <- rep("miRNA", length(sncInteract$regulator))
  } else {
    # do nothing
  }
  return(sncInteract)
}


#' Construct a transcriptional regulatory network
#'
#' @aliases constructTRN
#'
#' @description This function constructs a transcriptional regulatory network
#'              based on the omics data and the external interaction data
#'              provided. The network is stored as an igraph object in the
#'              TRNet S4 class, which is a mRNA-TF-smallRNA network.
#'
#' @param objMOList A MOList object containing all omics data that the user
#'                  wants to use to construct the network
#' @param omiCutoffs A OMICutoffs object containing the cutoffs for the omics
#'                   data
#' @param useOmics A character vector containing the omics data that the user
#'                 wants to use to construct the network. The available omics
#'                 data are "RNAseq", "smallRNAseq", "proteomics", "ATACseq",
#'                 and "extInteractions". By default, all the omics data will
#'                 be used.
#' @param smallRNAtypes A character vector containing the small RNA types that
#'                      the user wants to use to construct the network. The
#'                      available types are "miRNA", "piRNA", "snRNA", "snoRNA",
#'                      "tRNA", and "circRNA". If "all" is specified, then all
#'                      the small RNA types will be used.
#' @param targetDirection A single-length character vector indicating the
#'                        direction of expression change of the target genes.
#'                        Default to "both". Other options include "up" and
#'                        "down".
#' @param predicted A logical value indicating whether to perform predicted
#'                  inference of small RNA - mRNA interactions. Default to TRUE.
#'                  Prediction utilizes the correlation of expression between
#'                  small RNAs and mRNAs, which is highly recommended for adding
#'                  an additional level of condition-specificity to the network.
#' @param ntree The number of trees to grow in the random forest model, Must be
#'              a positive integer. Default is 1000.
#' @param nthreads The number of threads to use for parallel processing. Must be
#'                 a positive integer between 1 and the number of available
#'                 cores on the user's machine. Default is 1.
#' @param treeMethod The method to use for the random forest model. Either "RF"
#'                   or "ET" for random forest and extra trees, respectively.
#'                   See the documentation for the \code{\link[GENIE3]{GENIE3}}
#'                   package for details. Default is "RF".
#' @param seed The seed to use for the random forest model. Must be a positive
#'             integer. Default is 91711.
#'
#' @return A TRNet object containing the transcriptional regulatory network
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom Rdpack reprompt
#' @importFrom igraph graph_from_data_frame
#' @importFrom GENIE3 GENIE3
#'
#' @references
#' \insertRef{dplyr}{IntegraTRN}
#'
#' \insertRef{csardi2006igraph}{IntegraTRN}
#'
#' \insertRef{huynh2010inferring}{IntegraTRN}
#'
#' @examples
#' # Use the package data
#' data("expMOList")
#'
#' # Set the cutoffs for the omics data
#' omiCutoffs <- setOmicCutoffs() # use default cutoffs
#'
#' # Construct the network
#' \dontrun{
#' myTRNet <- constructTRN(expMOList, omiCutoffs, targetDirection = "up")
#' plotNetwork(myTRNet)
#' }
#'
constructTRN <- function(objMOList,
                         omiCutoffs,
                         useOmics = c(
                           "RNAseq", "smallRNAseq",
                           "proteomics", "ATACseq",
                           "extInteractions"
                         ),
                         smallRNAtypes = "all",
                         targetDirection = c("up", "down", "both"),
                         predicted = TRUE,
                         ntree = 1000,
                         nthreads = 1,
                         treeMethod = "RF",
                         seed = 91711) {
  # Check that the selected omics data are valid strings
  if (!all(useOmics %in% c(RNA, SMALLRNA, PROTEIN, ATAC, "extInteractions"))) {
    stop("Invalid omics data specification.")
  } else {
    # Continue
  }

  # Check the available omics data
  rnaSeq <- is.null(objMOList$DERNAseq) || (!"RNAseq" %in% useOmics)
  smallRNAseq <- is.null(objMOList$DEsmallRNAseq) ||
    (!"smallRNAseq" %in% useOmics)
  proteomics <- is.null(objMOList$DEproteomics) ||
    (!"proteomics" %in% useOmics)
  atacSeq <- is.null(objMOList$DEATAC) || (!"ATACseq" %in% useOmics)
  extTF2gene <- is.null(objMOList$extInteractions$tf2Genes) ||
    (!"extInteractions" %in% useOmics)
  extmiR2gene <- is.null(objMOList$extInteractions$miR2Genes) ||
    (!"extInteractions" %in% useOmics)
  checkOmics <- c(
    rnaSeq, smallRNAseq, proteomics, atacSeq, extTF2gene,
    extmiR2gene
  )

  # Validate inputs
  if (rnaSeq) {
    stop("No differential RNAseq data provided. Must perform differential
    analysis before constructing the network.")
  } else if (length(targetDirection) != 1 ||
    !targetDirection %in% c("up", "down", "both")) {
    stop("The targetDirection must be one of \"up\", \"down\" or \"both\".")
  } else if (!all(smallRNAtypes %in% c(SMALLRNA_CATEGORIES, "all"))) {
    stop("Invalid small RNA type specification. See ?constructTRN for details.")
  } else {
    # Do nothing
  }

  # Logic for validating whether a small RNA - TF - mRNA network can be
  # constructed
  if ((sum(checkOmics) == 5 && !rnaSeq) ||
    (sum(checkOmics) == 4 && !rnaSeq && !proteomics)) {
    # 1. Only RNAseq is available or only RNAseq and proteomics are available
    stop(paste0(
      "A small RNA - TF - mRNA network cannot be constructed with ",
      "only RNAseq and/or proteomics data. It is highly recommended to ",
      "search for externally curated interactions based on the RNAseq ",
      "target genes."
    ))
  } else if (extmiR2gene && extTF2gene && !predicted) {
    # 2. No external data and predictive inference is disabled
    stop(paste0(
      "A small RNA-TF-mRNA network cannot be constructed without ",
      "any external data and predictive inference disabled. Without any ",
      "externally curated global interaction data, the network only relies ",
      "on the predictive inference of regulator -target interactions."
    ))
  }

  # Based on the availability of the data, different methods will be used to
  # construct the network

  # Omic data used
  omics <- RNA

  # ====== PART 1: Filter target genes by proteomics data ======================
  # Only genes that shows consistent expression changes in both RNAseq and
  # proteomics data will be used to construct the network
  # Note: target gene direction is initially filtered at rnaTOPTag construction
  rnaTopTag <- TOPTag(objMOList$DERNAseq,
    logFCCutoff = omiCutoffs$rnaLogFC,
    pCutoff = omiCutoffs$rnaAdjPval,
    topGenes = omiCutoffs$rnaTopGenes,
    direction = targetDirection
  )
  if (proteomics || is.null(objMOList$gene2protein)) {
    warning(paste0(
      "No differential proteomics data or gene to protein ",
      "mapping. Skipping filtering target genes by proteomics data."
    ))
  } else {
    # Filter the target genes by proteomics data
    protTopTag <- TOPTag(objMOList$DEproteomics,
      logFCCutoff = omiCutoffs$proteomicsLogFC,
      pCutoff = omiCutoffs$proteomicsAdjPval,
      topGenes = 1,
      direction = targetDirection
    )
    # Filter the target genes by proteomics data
    rnaTopTag <- filterTargetGenes(
      rnaTopTag,
      protTopTag,
      objMOList$gene2protein
    )
    omics <- union(omics, PROTEIN)
  }

  # ====== PART 2: small RNA - mRNA interactions ===============================
  if (all(smallRNAtypes == "all")) {
    # Use all the small RNA types
    smallRNAtypes <- SMALLRNA_CATEGORIES
  } else if (all(smallRNAtypes %in% SMALLRNA_CATEGORIES)) {
    # Do nothing
  } else {
    stop("Invalid small RNA type specification.")
  }
  # check only if smallRNAseq is used
  if (!is.null(objMOList$annoSncRNA) &&
    length(objMOList$annoSncRNA) == 1 &&
    objMOList$annoSncRNA == HUMAN) {
    sncAnno <- SNCANNOLIST_HSAPIENS
  } else {
    sncAnno <- objMOList$annoSncRNA
  }

  # Perform predicted inference of small RNA - mRNA interactions, if specified
  if (predicted == TRUE) {
    if (smallRNAseq) {
      warning(paste0(
        "No differential smallRNAseq data provided. Skipping ",
        "predicted inference of small RNA - mRNA interactions."
      ))
      smallRNAmRNAPred <- NULL
    } else {
      smallRNAmRNAPred <- predictSncmRNAInteractions(
        objMOList = objMOList,
        omiCutoffs = omiCutoffs,
        rnaTopTag = rnaTopTag,
        smallRNATypes = smallRNAtypes,
        ntree = ntree,
        nthreads = nthreads,
        treeMethod = treeMethod,
        seed = seed
      )
      smallRNAmRNAPred <- lapply(smallRNAmRNAPred, function(x) as.character(x))
      omics <- union(omics, SMALLRNA)
    }
  } else {
    smallRNAmRNAPred <- NULL
  }
  # User imported small RNA - mRNA interactions
  if (extmiR2gene) {
    # No user-imported small RNA - mRNA interactions
    miRNAmRNAExt <- NULL
  } else {
    miRNAmRNAExt <- objMOList$extInteractions$miR2Genes
    omics <- union(omics, SMALLRNA)
  }

  if (smallRNAseq) {
    deSmallRNA <- NULL
  } else {
    deSmallRNA <- exportDE(objMOList$DEsmallRNAseq)
  }

  # Finalize small RNA - mRNA interactions
  smallRNAmRNA <- combineSncInteractions(
    predInteract = smallRNAmRNAPred,
    extInteract = miRNAmRNAExt,
    sncAnno = sncAnno,
    deResultRNA = objMOList$DERNAseq %>% exportDE(),
    deResultSmallRNA = deSmallRNA,
    smallRNATypes = smallRNAtypes
  )

  # ====== PART 3: TF - mRNA interactions ======================================
  # Check if the user has provided TF - mRNA interactions
  if (extTF2gene) {
    # No user-imported TF - mRNA interactions, skip TF - mRNA interactions
    tfmRNA <- NULL
  } else {
    extTFmRNA <- objMOList$extInteractions$tf2Genes
    extTFmRNA$regulatorType <- rep("TF", length(extTFmRNA$regulator))
    # If ATACseq data is available, filter TFs by enriched motifs
    if (atacSeq) {
      # Skip filtering silently since no ATACseq data
      tfmRNA <- extTFmRNA
    } else if (!inherits(objMOList$DEATAC, "PEAKTag")) {
      warning(paste0(
        "ATACseq data has no annotated motifs. See ?annotateATACPeaksMotif ",
        "for details. Skipping filtering TFs by ATACseq data."
      ))
      tfmRNA <- extTFmRNA
    } else {
      # Filter TFs by ATACseq data
      sel <- selectedMotifs(
        enrichedMotifs = objMOList$DEATAC$motifEnrichment,
        pValueAdj = omiCutoffs$atacMotifAdjPval,
        pValue = omiCutoffs$atacMotifPval,
        log2FEnrich = omiCutoffs$atacMotifLogFC
      )
      validTFs <- motifNames(objMOList$DEATAC$motifEnrichment[sel, ])
      # Filter to retain only valid interactions
      sel <- extTFmRNA$regulator %in% validTFs
      tfmRNA <- lapply(extTFmRNA, function(x) x[sel])
      omics <- union(omics, ATAC)
    }
  }

  # ====== PART 4: Combine all the interactions ================================
  # Combine all the interactions
  trnInteractions <- rbind(data.frame(smallRNAmRNA), data.frame(tfmRNA))
  # Access any potential duplicated interactions (potentially user-imported)
  duplicatedInt <- duplicated(trnInteractions[, NETWORK_FIELD[1:2]])
  trnInteractions <- trnInteractions[!duplicatedInt, ]
  # Filter by target
  trnInteractions <- trnInteractions %>%
    dplyr::filter(target %in% rownames(rnaTopTag %>% exportDE()))
  # Construct the network
  trnNetwork <- TRNet(
    TRNmetadata = trnInteractions,
    predicted = predicted,
    omics = omics
  )
  return(trnNetwork)
}


# [END]
