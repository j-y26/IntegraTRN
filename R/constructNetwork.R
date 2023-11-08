# Purpose: Constructing transcriptional regulatory network
# Author: Jielin Yang
# Date: 2023-11-07
# Version: 1.0
# Bugs and Issues: None


#' Loading external interaction data into the MOList object
#'
#' @description This function loads external interaction data into the MOList
#'              object. The external data can be curated non-tissue/sample
#'              specific interaction data, or any interactions that the users
#'              has determined to encompass a global picture of the regulatory
#'              interactions of the target genes.
#'
#' @param objMOList A MOList object containing the omics data
#' @param upregGenes2miR A data frame containing the upregulated genes and their
#'                       regulatory miRNAs. Must contain the columns "ID"
#'                       of the miRNAs and "Target" of the upregulated genes.
#' \itemize{
#'  \item \code{ID}: A column containing the miRNA names
#'  \item \code{Target}: A column containing the target genes
#' }
#' @param downregGenes2miR A data frame containing the downregulated genes and
#'                         their regulatory miRNAs. See above format.
#' @param upregGenes2TF A data frame containing the upregulated genes and their
#'                      regulatory transcription factors. See above format.
#' @param downregGenes2TF A data frame containing the downregulated genes and
#'                        their regulatory transcription factors. See above
#'                        format.
#'
#' @return An object of class MOList, with a element "extInteractions" added,
#'         which is a list of data frames
#' \itemize{
#' \item \code{upregGenes2miR}: A data frame containing the upregulated genes
#'                              and their regulatory miRNAs
#' \item \code{downregGenes2miR}: A data frame containing the downregulated
#'                                genes and their regulatory miRNAs
#' \item \code{upregGenes2TF}: A data frame containing the upregulated genes
#'                             and their regulatory transcription factors
#' \item \code{downregGenes2TF}: A data frame containing the downregulated
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
    stop("Please provide both upregulated and downregulated genes to miRNA
    interactions.")
  } else if (xor(is.null(upregGenes2TF), is.null(downregGenes2TF))) {
    stop("Please provide both upregulated and downregulated genes to TF
    interactions.")
  } else {
    # Continue
  }

  if (!is.null(upregGenes2miR)) {
    if (!(all(c("ID", "Target") %in% colnames(upregGenes2miR)) &&
      all(c("ID", "Target") %in% colnames(downregGenes2miR)))) {
      stop("Invalid column names provided. Must contain \"ID\" and \"Target\",
      see ?loadExtInteractions for details.")
    } else {
      # Do nothing
    }
  } else if (!is.null(upregGenes2TF)) {
    if (!(all(c("ID", "Target") %in% colnames(upregGenes2TF)) &&
      all(c("ID", "Target") %in% colnames(downregGenes2TF)))) {
      stop("Invalid column names provided. Must contain \"ID\" and \"Target\",
      see ?loadExtInteractions for details.")
    } else {
      # Do nothing
    }
  } else {
    # Do nothing
  }

  # Setting the external interactions to the MOList object
  objMOList$extInteractions <- list(
    upregGenes2miR = upregGenes2miR[, INTERACTION_FIELDS],
    downregGenes2miR = downregGenes2miR[, INTERACTION_FIELDS],
    upregGenes2TF = upregGenes2TF[, INTERACTION_FIELDS],
    downregGenes2TF = downregGenes2TF[, INTERACTION_FIELDS]
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
#'              TFs will be used to construct the network.
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
                           proteomicsLogFC = 1) {
  # Check the inputs
  inputs <- c(rnaAdjPval, rnaLogFC, rnaTopGenes,
              smallRNAAdjPval, smallRNALogFC, smallRNATopGenes,
              proteomicsAdjPval, proteomicsLogFC)
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
    proteomicsLogFC = proteomicsLogFC
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
}


#' Construct a transcriptional regulatory network
#' 
#' @description This function constructs a transcriptional regulatory network
#'              based on the omics data and the external interaction data
#'              provided. The network is stored as an igraph object in the
#'              TRNet S4 class.
#' 
#' @param objMOList A MOList object containing all omics data that the user
#'                  wants to use to construct the network
#' 
#' 
constructTRN <- function(objMOList
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
  
  # Based on the availability of the data, different methods will be used to
  # construct the network
  

}