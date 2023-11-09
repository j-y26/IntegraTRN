# Purpose: Constructing transcriptional regulatory network
# Author: Jielin Yang
# Date: 2023-11-07
# Version: 1.0
# Bugs and Issues: None


# Define global variables
NETOWRK_FIELD <- c("regulator", "target", "regulatorType")


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

  if (!all(NETOWRK_FIELD[1:2] %in% names(adjList))) {
    warning("Invalid element names provided.")
    warning("Setting the first element as regulator and second as target.")
    names(adjList)[1:2] <- NETOWRK_FIELD[1:2]
  } else {
    # Do nothing
  }
  return(adjList[, NETOWRK_FIELD])
}


#' Loading external interaction data into the MOList object
#'
#' @description This function loads external interaction data into the MOList
#'              object. The external data can be curated non-tissue/sample
#'              specific interaction data, or any interactions that the users
#'              has determined to encompass a global picture of the regulatory
#'              interactions of the target genes.
#'
#' @param objMOList A MOList object containing the omics data
#' @param upregGenes2miR A list containing the upregulated genes and their
#'                       regulatory miRNAs. Must contain elements "regulator"
#'                       for the miRNAs and "target" for the upregulated genes.
#' \itemize{
#'  \item \code{regulator}: A character vector containing the miRNAs
#'  \item \code{target}: A character vector containing the upregulated genes
#' }
#' @param downregGenes2miR A list containing the downregulated genes and
#'                         their regulatory miRNAs. See above format.
#' @param upregGenes2TF A list containing the upregulated genes and their
#'                      regulatory transcription factors. See above format.
#' @param downregGenes2TF A list containing the downregulated genes and
#'                        their regulatory transcription factors. See above
#'                        format.
#'
#' @return An object of class MOList, with a element "extInteractions" added,
#'         which is a list of lists. Each lower level list must follow the
#'         format of the input data, i.e., must contain elements "regulator"
#'        and "target". See the examples for details.
#' \itemize{
#' \item \code{upregGenes2miR}: A list containing the upregulated genes
#'                              and their regulatory miRNAs
#' \item \code{downregGenes2miR}: A list containing the downregulated
#'                                genes and their regulatory miRNAs
#' \item \code{upregGenes2TF}: A list containing the upregulated genes
#'                             and their regulatory transcription factors
#' \item \code{downregGenes2TF}: A list containing the downregulated
#'                               genes and their regulatory transcription
#'                               factors
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
#' upregGenes2miR <- data.frame(
#'   ID = paste0("miR_", seq_len(10)),
#'   Target = paste0("gene_", seq_len(10))
#' )
#' downregGenes2miR <- data.frame(
#'   ID = paste0("miR_", seq(11, 20)),
#'   Target = paste0("gene_", seq(11, 20))
#' )
#' upregGenes2TF <- data.frame(
#'   ID = paste0("TF_", seq(21, 30)),
#'   Target = paste0("gene_", seq(21, 30))
#' )
#' downregGenes2TF <- data.frame(
#'   ID = paste0("TF_", seq(31, 40)),
#'   Target = paste0("gene_", seq(31, 40))
#' )
#'
#' # Load the external interaction data into the MOList object
#' myMOList <- loadExtInteractions(
#'   myMOList,
#'   upregGenes2miR = upregGenes2miR,
#'   downregGenes2miR = downregGenes2miR,
#'   upregGenes2TF = upregGenes2TF,
#'   downregGenes2TF = downregGenes2TF
#' )
#'
loadExtInteractions <- function(objMOList,
                                upregGenes2miR = NULL,
                                downregGenes2miR = NULL,
                                upregGenes2TF = NULL,
                                downregGenes2TF = NULL) {
  # Check the input data
  if (is.null(upregGenes2miR) && is.null(downregGenes2miR) &&
    is.null(upregGenes2TF) && is.null(downregGenes2TF)) {
    stop("Please provide at least one type of the interaction data.")
  } else if (xor(is.null(upregGenes2miR), is.null(downregGenes2miR))) {
    # Both upregulated and downregulated genes must be provided, if any
    # miRNA interactions are provided
    stop("Please provide both upregulated and downregulated genes to miRNA
    interactions.")
  } else if (xor(is.null(upregGenes2TF), is.null(downregGenes2TF))) {
    # Same logic for TF interactions
    stop("Please provide both upregulated and downregulated genes to TF
    interactions.")
  } else {
    # Continue
  }

  if (!is.null(upregGenes2miR)) {
    # In case where the names of the elements of the list are not provided
    # correctly, set the first element to be the regulator and the second
    # element to be the target
    upregGenes2miR <- validateInteractionAdjList(upregGenes2miR)
    downregGenes2miR <- validateInteractionAdjList(downregGenes2miR)
  } else {
    # Do nothing
  }
  if (!is.null(upregGenes2TF)) {
    upregGenes2TF <- validateInteractionAdjList(upregGenes2TF)
    downregGenes2TF <- validateInteractionAdjList(downregGenes2TF)
  } else {
    # Do nothing
  }

  # Setting the external interactions to the MOList object
  objMOList$extInteractions <- list(
    upregGenes2miR = upregGenes2miR,
    downregGenes2miR = downregGenes2miR,
    upregGenes2TF = upregGenes2TF,
    downregGenes2TF = downregGenes2TF
  )

  return(objMOList)
}


#' Setting cutoffs for the omics data
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
#' @param rnaAdjPval The adjusted p-value cutoff for the RNAseq data
#' @param rnaLogFC The log fold change cutoff for the RNAseq data
#' @param rnaTopGenes A numeric value indicating either the fraction (0-1) of 
#'                    top differential genes or the number (1-Inf) of top
#'                    differential genes. If the number specified is greater
#'                    than the number of DE genes based on logFCCutoff and
#'                    pCutoff, then topGenes will be set to the number of DE
#'                    genes. This is to select the top differential mRNAs.
#'                    To include all the differential mRNAs, set this to 1.
#' @param smallRNAAdjPval The adjusted p-value cutoff for the smallRNAseq data
#' @param smallRNALogFC The log fold change cutoff for the smallRNAseq data
#' @param smallRNATopGenes A numeric value indicating either the fraction (0-1)
#'                         of top differential genes or the number (1-Inf) of
#'                         top differential genes. If the number specified is
#'                         greater than the number of DE genes based on
#'                         logFCCutoff and pCutoff, then topGenes will be set
#'                         to the number of DE genes. This is to select the top
#'                         differential small RNAs. To include all the 
#'                         differential small RNAs, set this to 1.
#' @param proteomicsAdjPval The adjusted p-value cutoff for the proteomics data
#' @param proteomicsLogFC The log fold change cutoff for the proteomics data
#' @param atacMotifAdjPval The adjusted p-value cutoff for the ATACseq motif
#'                         enrichment analysis
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
#' omiCutoffs <- setOmicCutoffs(rnaAdjPval = 0.05, 
#'                              rnaLogFC = 2,
#'                              rnaTopGenes = 0.1,      # Top 10% of DE genes
#'                              smallRNAAdjPval = 0.05,
#'                              smallRNALogFC = 2,
#'                              smallRNATopGenes = 200, # Top 200 DE small RNAs
#'                              proteomicsAdjPval = 0.05,
#'                              proteomicsLogFC = 1)
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
                           rnaLogFC = 1,
                           rnaTopGenes = 1,
                           smallRNAAdjPval = 0.05,
                           smallRNALogFC = 1,
                           smallRNATopGenes = 1,
                           proteomicsAdjPval = 0.05,
                           proteomicsLogFC = 1,
                           atacMotifAdjPval = 0.05,
                           atacMotifPval = NULL,
                           atacMotifLogFC = NULL) {
  # Check the inputs
  inputs <- c(rnaAdjPval, rnaLogFC, rnaTopGenes,
              smallRNAAdjPval, smallRNALogFC, smallRNATopGenes,
              proteomicsAdjPval, proteomicsLogFC,
              atacMotifAdjPval)
  if (!all(is.numeric(inputs))) {
    stop("The cutoffs must be numeric.")
  } else if (!all(inputs >= 0)) {
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
#' @description This function prints the cutoffs for the omics data.
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
#' 
#' # The output will be:
#' # RNAseq adjusted p-value cutoff: 0.05
#' # RNAseq log fold change cutoff: 1
#' # RNAseq selecting top 100% of DE genes
#' # smallRNAseq adjusted p-value cutoff: 0.05
#' # smallRNAseq log fold change cutoff: 1
#' # smallRNAseq selecting top 100% of DE genes
#' # Proteomics adjusted p-value cutoff: 0.05
#' # Proteomics log fold change cutoff: 1
#' 
print.OMICutoffs <- function(x, ...) {
  cat("RNAseq adjusted p-value cutoff:", x$rnaAdjPval, "\n")
  cat("RNAseq log fold change cutoff:", x$rnaLogFC, "\n")
  if (x$rnaTopGenes <= 1) {
    cat("RNAseq selecting top", x$rnaTopGenes * 100, "% of DE genes", "\n")
  } else {
    cat("RNAseq selecting top", x$rnaTopGenes, "DE genes", "\n")
  }
  cat("smallRNAseq adjusted p-value cutoff:", x$smallRNAAdjPval, "\n")
  cat("smallRNAseq log fold change cutoff:", x$smallRNALogFC, "\n")
  if (x$smallRNATopGenes <= 1) {
    cat("smallRNAseq selecting top", x$smallRNATopGenes * 100, "% of DE genes", 
    "\n")
  } else {
    cat("smallRNAseq selecting top", x$smallRNATopGenes, "DE genes", "\n")
  }
  cat("Proteomics adjusted p-value cutoff:", x$proteomicsAdjPval, "\n")
  cat("Proteomics log fold change cutoff:", x$proteomicsLogFC, "\n")
  if (is.null(x$atacMotifPval)) {
    cat("ATACseq motif adjusted p-value cutoff:", x$atacMotifAdjPval, "\n")
  } else {
    cat("ATACseq motif p-value cutoff:", x$atacMotifPval, "\n")
  }
  cat("ATACseq motif log fold enrichment cutoff:", x$atacMotifLogFC, "\n")
}




#' Select miRNA inverse correlation
#' 
#' @description This function selects the miRNAs in which its expression is
#'              inversely correlated with the expression of the target gene.
#' 
#' @param exprAdjList A list representing the regulator-target interactions
#' \itemize{
#' \item \code{regulator}: A character vector containing the regulators
#' \item \code{target}: A character vector containing the targets
#' }
#' @param rnaDETag A DETag object containing the differential RNAseq data
#' @param smallDETag A DETag object containing the differential smallRNAseq data
#' @param smallRNATypes A character vector containing the small RNA types that
#'                      the user wants to use to construct the network. The
#'                      available types are "miRNA", "piRNA", "snRNA", "snoRNA",
#'                      "tRNA", and "circRNA".
#' 
#' @return A list with filtered regulator-target interactions
#' 
filtermiRNAinverseCorr <- function(exprAdjList, 
                                   deResultRNA, 
                                   deResultSmallRNA,
                                   smallRNATypes = SMALLRNA_CATEGORIES) {
  if ("miRNA" %in% smallRNATypes) {
      # Extract only inverse small RNA - mRNA interactions
      # i.e., downregulated small RNA - upregulated mRNA
      genesmiRNA <- extractDirectionalGenes(smallDETag) %>%
        intersect(., sncAnno[["miRNA"]])
      genesRNA <- extractDirectionalGenes(rnaDETag)

      downregmiRNA <- genesmiRNA$down
      upregmiRNA <- genesmiRNA$up
      downregRNA <- genesRNA$down
      upregRNA <- genesRNA$up

      # Extract the predicted interactions
      inverseRelation <- (exprAdjList$regulator %in% downregmiRNA &
                          exprAdjList$target %in% upregRNA) |
                        (exprAdjList$regulator %in% upregmiRNA &
                          exprAdjList$target %in% downregRNA)
      notmiRNA <- !(exprAdjList$regulator %in% sncAnno[["miRNA"]])
      exprAdjList <- lapply(exprAdjList, 
                              function(x) x[inverseRelation | notmiRNA])
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
#' @rnaTopTag A TOPTag object containing the top differential genes from the
#'            RNAseq data
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
  smallRNATopTag <- TOPTag(smallDETag,
                           logFCCutoff = omiCutoffs$smallRNALogFC,
                           pCutoff = omiCutoffs$smallRNAAdjPval,
                           topGenes = omiCutoffs$smallRNATopGenes,
                           direction = "both")

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
  intersectedList <- lapply(adjList1, function(x) x[x$pair %in% commonInteract])
  intersectedList$pair <- NULL

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
filterTargetGenes <- function(rnaTopTag, protTopTag, mapping) {
  # Filter the target genes by proteomics data
  validProteins <- protTopTag %>% exportDE() %>% rownames()
  validTargets <- mapping$gene[mapping$protein %in% validProteins]

  # Redefine the target genes
  rnaTopTag <- filterGenes(rnaTopTag, validTargets)
}




#' Construct a transcriptional regulatory network
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
#' @param smallRNATypes A character vector containing the small RNA types that
#'                      the user wants to use to construct the network. The
#'                      available types are "miRNA", "piRNA", "snRNA", "snoRNA",
#'                     "tRNA", and "circRNA". If "all" is specified, then all
#'                     the small RNA types will be used.
#' @param targetDirection A character vector indicating the direction of the
#'                        target genes. Default to "both". Other options include
#'                       "up" and "down".
#' @param predicted A logical value indicating whether to perform predicted
#'                  inference of small RNA - mRNA interactions. Default to TRUE.
#' @param ntree The number of trees to grow in the random forest model
#' @param nthreads The number of threads to use for parallel processing
#' @param treeMethod The method to use for the random forest model. Either "RF"
#'                   or "ET". See the documentation for the GENIE3 package for
#'                   details.
#' @param seed The seed to use for the random forest model
#' 
#' 
constructTRN <- function(objMOList, 
                         omiCutoffs,
                         smallRNAtypes = "all",
                         targetDirection = c("up", "down", "both"),
                         predicted = TRUE,
                         ntree = 1000,
                         nthreads = 1,
                         treeMethod = "RF",
                         seed = 91711
                         ) {
  # Check the available omics data
  rnaSeq <- is.null(objMOList$DERNAseq)
  smallRNAseq <- is.null(objMOList$DEsmallRNAseq)
  proteomics <- is.null(objMOList$DEProteomics)
  atacSeq <- is.null(objMOList$DEATAC)
  extTF2gene <- is.null(objMOList$extInteractions$upregGenes2TF) &&
    is.null(objMOList$extInteractions$downregGenes2TF)
  extmiR2gene <- is.null(objMOList$extInteractions$upregGenes2miR) &&
    is.null(objMOList$extInteractions$downregGenes2miR)

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

  # Based on the availability of the data, different methods will be used to
  # construct the network

  # Omic data used
  omics <- RNA

  # ====== PART 1: Filter target genes by proteomics data ======================
  # Only genes that shows consistent expression changes in both RNAseq and
  # proteomics data will be used to construct the network
  if (proteomics || is.null(objMOList$gene2protein)) {
    warning("No differential proteomics data or gene to protein mapping.")
    warning("Skipping filtering target genes by proteomics data.")
    rnaTopTag <- TOPTag(objMOList$DERNAseq,
                        logFCCutoff = omiCutoffs$rnaLogFC,
                        pCutoff = omiCutoffs$rnaAdjPval,
                        topGenes = omiCutoffs$rnaTopGenes,
                        direction = targetDirection)
  } else {
    # Filter the target genes by proteomics data
    rnaTopTag <- TOPTag(objMOList$DERNAseq,
                        logFCCutoff = omiCutoffs$rnaLogFC,
                        pCutoff = omiCutoffs$rnaAdjPval,
                        topGenes = omiCutoffs$rnaTopGenes,
                        direction = targetDirection)
    protTopTag <- TOPTag(objMOList$DEProteomics,
                         logFCCutoff = omiCutoffs$proteomicsLogFC,
                         pCutoff = omiCutoffs$proteomicsAdjPval,
                         topGenes = 1,
                         direction = targetDirection)
    # Filter the target genes by proteomics data
    rnaTopTag <- filterTargetGenes(rnaTopTag, protTopTag)
    omic <- union(omic, PROTEIN)
  }

  # ====== PART 2: small RNA - mRNA interactions ===============================
  if (smallRNAtypes == "all") {
    # Use all the small RNA types
    smallRNAtypes <- SMALLRNA_CATEGORIES
  } else {
    # Do nothing
  }

  # Perform predicted inference of small RNA - mRNA interactions, if specified
  if (predicted == TRUE) {
    if (smallRNAseq) {
      warning("No differential smallRNAseq data provided. Skipping predicted
      inference of small RNA - mRNA interactions.")
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
      omic <- union(omic, SMALLRNA)
    }
  } else {
    smallRNAmRNAPred <- NULL
  }
  # User imported small RNA - mRNA interactions
  if (extmiR2gene) {
    # No user-imported small RNA - mRNA interactions
    smallRNAmRNAExt <- NULL
  } else {
    smallRNAmRNAExt <- switch(targetDirection,
      "up" = objMOList$extInteractions$upregGenes2miR,
      "down" = objMOList$extInteractions$downregGenes2miR,
      "both" = mapply(c, objMOList$extInteractions$upregGenes2miR,
                      objMOList$extInteractions$downregGenes2miR,
                      SIMPLIFY = FALSE)
    )
    omic <- union(omic, SMALLRNA)
  }

  # Finalize small RNA - mRNA interactions
  if (is.null(smallRNAmRNAPred) && is.null(smallRNAmRNAExt)) {
    # No small RNA - mRNA interactions
    smallRNAmRNA <- NULL
  } else {
    # Intersect the predicted and user-imported interactions
    smallRNAmRNA <- intersectInteractions(smallRNAmRNAPred, smallRNAmRNAExt)
    # Filter for inverse correlation on miRNA - mRNA interactions
    smallRNAmRNA <- filtermiRNAinverseCorr(smallRNAmRNA,
                                           objMOList$DERNAseq,
                                           objMOList$DEsmallRNAseq,
                                           smallRNAtypes)
    # Annotate the small RNA types
    if (all(objMOList$annoSncRNA == HUMAN)) {
      sncAnno <- SNCANNOLIST_HSAPIENS
    } else {
      sncAnno <- objMOList$annoSncRNA
    }
    smallRNAmRNA$regulatorType <- sapply(
      smallRNAmRNA$regulator,
      findGeneType, 
      annotation = sncAnno
    )
  }

  # ====== PART 3: TF - mRNA interactions ======================================
  # Check if the user has provided TF - mRNA interactions
  if (extTF2gene) {
    # No user-imported TF - mRNA interactions, skip TF - mRNA interactions
    transcriptionFactormRNA <- NULL
  } else {
    extTFmRNA <- switch(targetDirection,
      "up" = objMOList$extInteractions$upregGenes2TF,
      "down" = objMOList$extInteractions$downregGenes2TF,
      "both" = mapply(c, objMOList$extInteractions$upregGenes2TF,
                      objMOList$extInteractions$downregGenes2TF,
                      SIMPLIFY = FALSE)
    )
    # Filter by ATACseq-validated TFs if available
    if (atacSeq) {
      # Skip filtering TFs by ATACseq data silently
    } else if (!inherits(objMOList$DEATAC, "PEAKTag")) {
      warning("ATACseq data has no annotated motifs. See ?annotateATACPeaksMotif
      for details.")
      warning("Skipping filtering TFs by ATACseq data.")
    } else {
      # Filter TFs by ATACseq data
      
      
    }

  }

}