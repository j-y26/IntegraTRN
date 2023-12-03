# Purpose: Plotting differentially expressed omics data
# Author: Jielin Yang
# Date: 2023-11-02
# Version: 1.0
# Bugs and Issues: None


# Define the global variables
NOANNO <- "Unannotated"


#' Annotate by expression
#'
#' @keywords internal
#'
#' @description This function annotates whether a gene is up-regulated,
#'              down-regulated, or not differentially expressed
#'
#' @param deg A data frame containing the differential expression analysis
#'            results, must follow the described format in DETag-class
#' @param log2FC The cutoff for log2 fold change, a positive number
#' @param adjP The cutoff for adjusted p-value
#'
#' @return A data frame containing the differential expression analysis with
#'         the expression annotated in the "expr" column
#'
#' @importFrom dplyr mutate case_when %>%
#'
annoExpr <- function(deg, log2FC, adjP) {
  deg <- deg %>%
    dplyr::mutate(expr = dplyr::case_when(
      logFC >= log2FC & padj < adjP ~ "Up-regulated",
      logFC <= -log2FC & padj < adjP ~ "Down-regulated",
      TRUE ~ "Not DE"
    ))
  return(deg)
}


#' Generate a base volcano plot
#'
#' @keywords internal
#'
#' @description This function generates a base volcano plot for visualizing
#'             differentially expressed omics data based on log2 fold change
#'             and adjusted p-value.
#'
#' @param deg A data frame containing the differential expression analysis
#'            results, must follow the described format in DETag-class
#' @param log2FC The cutoff for log2 fold change.
#' @param adjP The cutoff for adjusted p-value.
#' @param title The title for the plot.
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 aes scale_color_manual guides guide_legend theme_bw xlab
#' @importFrom ggplot2 ylab theme element_blank element_text element_line
#' @importFrom ggplot2 ggplot geom_point ggtitle geom_text
#' @importFrom dplyr mutate case_when
#'
#' @references
#' \insertRef{villanueva2019ggplot2}{IntegraTRN}
#'
#' \insertRef{dplyr}{IntegraTRN}
#'
plotBaseVolcano <- function(deg,
                        log2FC = 1,
                        adjP = 0.05,
                        title = NULL) {
  # Generate the plot
  vPlot <- ggplot2::ggplot(deg, ggplot2::aes(x = logFC, y = -log(padj, 10))) +
    ggplot2::theme_classic() +
    ggplot2::xlab(expression("log"[2] * "FC")) +
    ggplot2::ylab(expression("-log"[10] * "adj.p-value")) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::geom_hline(yintercept = -log(adjP, 10), linetype = "dashed") +
    ggplot2::geom_vline(xintercept = c(-log2FC, log2FC), linetype = "dashed") +
    ggplot2::theme(
      axis.line = ggplot2::element_line(size = 0.7),
      axis.ticks = ggplot2::element_line(size = 0.4),
      axis.text = ggplot2::element_text(size = 10)
    )

  # Add the title
  if (!is.null(title)) {
    vPlot <- vPlot + ggplot2::ggtitle(title)
  }

  return(vPlot)
}


#' Volcano plot for visualizing differentially expressed mRNA
#'
#' @aliases plotVolcanoRNA
#'
#' @description This function generates a volcano plot for visualizing
#'              differentially expressed mRNA based on log2 fold change and
#'              adjusted p-value. The user can specify the cutoffs for log2 fold
#'              change and adjusted p-value. The user can also specify a list of
#'              genes to highlight in the plot. The colors for up- and down-
#'              regulated genes can be specified by the user, but a default
#'              color scheme is provided.
#'
#' @param objMOList An MOList object containing the differential expression
#'                  results for mRNA.
#' @param omic The omic type, one of "RNAseq", "smallRNAseq", or "proteomics"
#' @param log2FC The cutoff for log2 fold change. Default is 0.
#' @param adjP The cutoff for adjusted p-value. Default is 0.05.
#' @param upColor The color for up-regulated genes. Default is "firebrick3".
#' @param downColor The color for down-regulated genes. Default is "dodgerblue3"
#' @param title The title for the plot. Default is NULL.
#' @param highlight A vector of gene names to highlight in the plot. Default
#'                  is NULL.
#' @param pointSize The size of the points in the plot. Default is 1.
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 aes scale_color_manual guides guide_legend theme_bw xlab
#' @importFrom ggplot2 ylab theme element_blank element_text element_line
#' @importFrom ggplot2 ggplot geom_point ggtitle geom_text
#' @importFrom dplyr mutate case_when
#'
#' @references
#' \insertRef{villanueva2019ggplot2}{IntegraTRN}
#'
#' \insertRef{dplyr}{IntegraTRN}
#'
#' @export
#'
#' @examples
#' # Use the package-provided example data
#' data(expMOList)
#'
#' # Example 1: Generate the volcano plot by default parameters
#' plotVolcanoRNA(expMOList)
#'
#' # Example 2: Generate the volcano plot with custom parameters
#' plotVolcano(expMOList,
#'   omic = "RNAseq",
#'   log2FC = 0,
#'   adjP = 0.01,results
#'   upColor = "purple",
#'   downColor = "green",
#'   title = "Volcano plot for DE mRNA"
#' )
#'
plotVolcano <- function(objMOList,
                        omic = c("RNAseq", "smallRNAseq", "proteomics"),
                        log2FC = 0,
                        adjP = 0.05,
                        upColor = "firebrick3",
                        downColor = "dodgerblue3",
                        title = NULL,
                        highlight = NULL,
                        pointSize = 1) {
  # Retrieve the differential expression results
  omic <- match.arg(omic, c(RNA, SMALLRNA, PROTEIN))
  deResult <- switch(omic,
    "RNAseq" = objMOList$DERNAseq,
    "smallRNAseq" = objMOList$DEsmallRNAseq,
    "proteomics" = objMOList$DEproteomics,
    stop("Invalid omic type. Please specify one of 'RMAseq', 'smallRNAseq',
    or 'proteomics'.")
  )

  # Check if the differential expression results are available
  if (is.null(deResult)) {
    stop("No differential expression results for ", omic, ". Please run
  diffOmics() first. See ?diffOmics for details.")
  } else {
    # Continue
  }

  # Retrieve the differential expression results for mRNA
  deResult <- deResult %>% exportDE()

  # Annotate expression
  deResult <- annoExpr(deResult, log2FC, adjP)
  deResult$gene <- rownames(deResult)

  # Generate the base volcano plot
  vPlot <- plotBaseVolcano(deResult, log2FC, adjP, title)

  # Color the up- and down-regulated genes
  vPlot <- vPlot +
    ggplot2::geom_point(ggplot2::aes(color = expr), size = pointSize) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes =
        list(size = 2.5)
    )) +
    ggplot2::scale_color_manual(values = c(
      "Up-regulated" = upColor,
      "Down-regulated" = downColor,
      "Not DE" = "grey50"
    ))

  # Highlight the genes of interest by labeling their gene names
  if (!is.null(highlight)) {
    vPlot <- vPlot +
      ggplot2::geom_label(data = deResult %>% filter(gene %in% highlight),
                      aes(label = gene),
                      size = 3)
  } else {
    # Continue
  }

  # Annotate the numbers of up- and down-regulated genes
  upNum <- sum(deResult$logFC >= log2FC & deResult$padj < adjP)
  downNum <- sum(deResult$logFC <= -log2FC & deResult$padj < adjP)
  vPlot <- vPlot + ggplot2::geom_text(
    ggplot2::aes(x = max(deResult$logFC), y = -max(log(deResult$padj, 10))),
    label = upNum,
    color = upColor,
    hjust = 1,
    vjust = 1
  ) + ggplot2::geom_text(
    ggplot2::aes(x = min(deResult$logFC), y = -max(log(deResult$padj, 10))),
    label = downNum,
    color = downColor,
    hjust = 0,
    vjust = 1
  )

  return(vPlot)
}


#' Annotate the type of small RNA in the differential expression results
#'
#' @keywords internal
#'
#' @description This function annotates the type of small RNA in the
#'              differential expression results.
#'
#' @param deg A data frame containing the differential expression analysis.
#'            Format consistent with the output of exportDE()
#' @param annoList A list of annotations for small RNA.
#'
#' @return A data frame containing the differential expression analysis with
#'         the type of small RNA annotated.
#'
#' @importFrom dplyr mutate case_when %>%
#'
#' @references
#' \insertRef{dplyr}{IntegraTRN}
#'
annoSncList <- function(deg, annoList) {
  deg <- deg %>% dplyr::mutate(type = dplyr::case_when(
    rownames(deg) %in% annoList$miRNA ~ "miRNA",
    rownames(deg) %in% annoList$piRNA ~ "piRNA",
    rownames(deg) %in% annoList$snoRNA ~ "snoRNA",
    rownames(deg) %in% annoList$snRNA ~ "snRNA",
    rownames(deg) %in% annoList$tRNA ~ "tRNA",
    rownames(deg) %in% annoList$circRNA ~ "circRNA",
    TRUE ~ NOANNO
  ))

  return(deg)
}


#' Volcano plot for visualizing differentially expressed small RNAs
#'
#' @aliases plotVolcanoSmallRNA
#'
#' @description This function generates a volcano plot for visualizing
#'              differentially expressed small RNAs based on log2 fold change
#'              and adjusted p-value. The user can specify the cutoffs for log2
#'              fold change and adjusted p-value. The user can also specify a
#'              l ist of genes to highlight in the plot. The users can define
#'              a color scheme used for color-coding each type of small RNAs.
#'              The default color scheme is provided.
#'
#' @param objMOList An MOList object containing the differential expression
#'                  results for small RNAs.
#' @param log2FC The cutoff for log2 fold change. Default is 0.
#' @param adjP The cutoff for adjusted p-value. Default is 0.05.
#' @param colScheme A RColorBrewer color scheme for color-coding each type of
#'                  small RNAs. Default is "BuPu". See ?RColorBrewer::brewer.pal
#'                  for details.
#' @param title A character vector specifying the title for the plot. Default
#'              is NULL.
#' @param highlight A vector of gene names to highlight in the plot.
#'                  Default is NULL.
#' @param pointSize The size of the points in the plot. Default is 1.
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 aes scale_color_manual guides guide_legend theme_bw xlab
#' @importFrom ggplot2 ylab theme element_blank
#' @importFrom ggplot2 ggplot geom_point ggtitle geom_text
#' @importFrom dplyr mutate case_when
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
#'
#' @references
#' \insertRef{villanueva2019ggplot2}{IntegraTRN}
#'
#' \insertRef{dplyr}{IntegraTRN}
#'
#' \insertRef{RColorBrewer}{IntegraTRN}
#'
#' @examples
#' # Use the package-provided example data
#' data(expMOList)
#'
#' # Example 1: Generate the volcano plot by default parameters
#' plotVolcanoSmallRNA(expMOList)
#'
#' # Example 2: Generate the volcano plot with custom parameters
#' plotVolcanoSmallRNA(expMOList,
#'   log2FC = 0,
#'   adjP = 0.01,
#'   colScheme = "YlOrRd",
#'   title = "Volcano plot for DE small RNA"
#' )
#'
plotVolcanoSmallRNA <- function(objMOList,
                                log2FC = 0,
                                adjP = 0.05,
                                colScheme = "BuPu",
                                title = NULL,
                                highlight = NULL,
                                pointSize = 1) {
  if (is.null(objMOList$DEsmallRNAseq)) {
    stop("No differential expression results for small RNA. Please run
  diffOmics() first. See ?diffOmics for details.")
  } else if (is.null(objMOList$annoSncRNA)) {
    stop("No annotations for small RNA. Please run annotateSmallRNA() first.")
  } else {
    # Continue
  }

  # Retrieve the differential expression results for small RNA
  degSmallRNAseq <- objMOList$DEsmallRNAseq %>% exportDE()

  # Annotate the type of small RNA
  sncAnno <- objMOList$annoSncRNA
  if (all(sncAnno == HUMAN)) {
    degSmallRNAseq <- annoSncList(degSmallRNAseq,
      annoList = SNCANNOLIST_HSAPIENS
    )
  } else if (is.list(sncAnno)) {
    degSmallRNAseq <- annoSncList(degSmallRNAseq, sncAnno)
  } else {
    stop("Incompatible annotation for small RNA. Please check the annotation.")
  }
  degSmallRNAseq <- annoExpr(degSmallRNAseq, log2FC, adjP)

  # Type of expressed small RNAs:
  sncTypes <- unique(degSmallRNAseq$type)

  # Re-label non-DE small RNAs
  degSmallRNAseq$type[degSmallRNAseq$expr == "Not DE"] <- "Not DE"
  degSmallRNAseq$gene <- rownames(degSmallRNAseq)

  # Generate the base volcano plot
  vPlot <- plotBaseVolcano(degSmallRNAseq, log2FC, adjP, title)

  # Generate a set of colors for each type of small RNA
  sncColors <- RColorBrewer::brewer.pal(length(sncTypes), colScheme)
  sncColorList <- stats::setNames(sncColors, sncTypes)

  # Color each type of small RNA
  vPlot <- vPlot +
    ggplot2::geom_point(ggplot2::aes(color = type), size = pointSize) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes =
        list(size = 2.5)
    )) +
    ggplot2::scale_color_manual(values = c(
      sncColorList,
      "Not DE" = "grey50"
    ))

  # Highlight the genes of interest by labeling their gene names
  if (!is.null(highlight)) {
    vPlot <- vPlot +
      ggplot2::geom_label(data = degSmallRNAseq %>% filter(gene %in% highlight),
                      aes(label = gene),
                      size = 3)
  } else {
    # Continue
  }

  return(vPlot)
}


#' Principle Component Analysis (PCA) for count-based omics data
#'
#' @aliases countPCA
#'
#' @description This function performs PCA analysis for count-based omics data.
#'              The user can specify the title for the plot.
#'
#' @details Raw counts are normalized by RLE method and then transformed by
#'         variance stabilizing transformation (VST) before PCA analysis.
#'
#' @param matrix A numeric matrix containing the count-based omics data, must be
#'               raw count expression matrix, not normalized by any method
#' @param groupBy A vector of factors for grouping samples
#' @param batch A vector of factors for batch effect, default is NULL
#' @param title The title for the plot, default is NULL
#' @param col1 The color for the first group, default is "#440154"
#' @param col2 The color for the second group, default is "#fde725". If
#'             continuous grouping variable is used, the color gradient will be
#'             used instead of the two colors
#'
#' @importFrom ggplot2 aes scale_color_manual guides guide_legend theme_bw xlab
#' @importFrom ggplot2 ylab theme element_blank element_text
#' @importFrom ggplot2 ggplot geom_point ggtitle geom_text
#' @importFrom ggplot2 coord_fixed scale_color_gradient
#' @importFrom dplyr mutate case_when
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq
#' @importFrom DESeq2 varianceStabilizingTransformation plotPCA
#'
#' @return A ggplot object
#'
#' @references
#' \insertRef{villanueva2019ggplot2}{IntegraTRN}
#'
#' \insertRef{dplyr}{IntegraTRN}
#'
#' \insertRef{love2014moderated}{IntegraTRN}
#'
#' @export
#'
#' @examples
#' # Use the package-provided example data
#' data(expMOList)
#' data(RNAseq_heart_samples)
#'
#' # Extract the information needed for PCA analysis for example RNAseq data
#' matrixRNAseq <- getRawData(expMOList, "RNAseq")
#' groupBy <- RNAseq_heart_samples$Age
#'
#' # Example 1: Perform PCA analysis without batch effect correction
#' pcaPlot <- countPCA(matrixRNAseq, groupBy)
#'
#' pcaPlot
#'
#' # Example 2: Perform PCA analysis with batch effect correction
#' # Create an example batch effect
#' batch <- rep(c("Batch1", "Batch2"), each = 30)[1:53]
#' pcaPlot <- countPCA(matrixRNAseq, groupBy, batch = batch)
#'
#' pcaPlot
#'
countPCA <- function(matrix,
                     groupBy,
                     batch = NULL,
                     title = NULL,
                     col1 = "#440154",
                     col2 = "#fde725") {
  # Check if the input is a numeric matrix
  if (!is.matrix(matrix) || !is.numeric(matrix)) {
    stop("The input must be a numeric matrix.")
  } else {
    # Continue
  }

  # Raw count expression normalization and variance stabilizing transformation
  designList <- DESeqDesign(groupBy, batch)
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = matrix,
    colData = designList$colData,
    design = designList$design
  )
  dds <- DESeq2::DESeq(dds)
  vts <- DESeq2::varianceStabilizingTransformation(dds)

  # Perform PCA analysis
  pcaData <- DESeq2::plotPCA(vts, intgroup = "group", returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  # Further customization of the PCA plot using ggplot2
  pcaPlot <- ggplot2::ggplot(
    pcaData,
    ggplot2::aes(x = -PC1, y = PC2, color = group)
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::coord_fixed(ratio = 0.85) +
    ggplot2::theme_bw() +
    ggplot2::xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ggplot2::ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(
        color = "black",
        size = 11,
        hjust = 0.5
      ),
      legend.title = ggplot2::element_blank()
    )
  # Color sample points based on the grouping variable
  if (is.numeric(groupBy)) {
    pcaPlot <- pcaPlot + ggplot2::scale_color_gradient(
      low = col1,
      high = col2
    )
  } else {
    pcaPlot <- pcaPlot + ggplot2::scale_color_manual(values = c(col1, col2))
  }

  # Add the title
  if (!is.null(title)) {
    pcaPlot <- pcaPlot + ggplot2::ggtitle(title)
  }

  return(pcaPlot)
}


#' Plotting results for Principle Component Analysis (PCA) for each type of
#' small RNA
#'
#' @aliases plotSmallRNAPCAs
#'
#' @description This function generates a PCA plot for each type of small RNA
#'              based on the normalized expression matrix. This function is
#'              used to explore the contribution of each type of small RNA to
#'              separating samples based on their phenotypes.
#'
#' @param objMOList An MOList object containing the raw expression matrix for
#'                  small RNAs
#' @param batch A vector of factors for batch effect, default is NULL. Used to
#'              correct batch effect in the PCA analysis for count normalization
#' @param col1 The color for the first group, default is "#440154"
#' @param col2 The color for the second group, default is "#fde725". If
#'             continuous grouping variable is used, the color gradient will be
#'             used instead of the two colors
#'
#' @return A list of ggplot objects
#'
#' @importFrom ggplot2 aes scale_color_manual guides guide_legend theme_bw xlab
#' @importFrom ggplot2 ylab theme element_blank element_text
#' @importFrom ggplot2 ggplot geom_point ggtitle geom_text
#' @importFrom ggplot2 coord_fixed scale_color_gradient
#' @importFrom dplyr mutate case_when
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq
#' @importFrom DESeq2 varianceStabilizingTransformation plotPCA
#'
#' @export
#'
#' @references
#' \insertRef{villanueva2019ggplot2}{IntegraTRN}
#'
#' \insertRef{dplyr}{IntegraTRN}
#'
#' \insertRef{love2014moderated}{IntegraTRN}
#'
#' @examples
#' \dontrun{
#' # Loading example smallRNAseq data
#' data("expMOList")
#' data("smallRNAseq_heart") # small RNAseq count matrix
#' data("smallRNAseq_heart_samples") # small RNAseq sample information
#'
#' # Adding the smallRNAseq data to the example MOList object
#' expMOList <- MOList(expMOList,
#'                     smallRNAseq = smallRNAseq_heart,
#'                     smallRNAGroupBy = smallRNAseq_heart_samples$Age)
#'
#' # Annotate the small RNAs using the package-provided human annotation
#' expMOList <- annotateSmallRNA(expMOList, "human")
#'
#' # Generate a list of PCA plots for each type of small RNA
#' pcaPlotList <- plotSmallRNAPCAs(expMOList)
#'
#' # Plot the PCA plot for miRNA
#' pcaPlotList$miRNA
#'
#' # Plot the PCA plot for piRNA
#' pcaPlotList$piRNA
#' }
#'
plotSmallRNAPCAs <- function(objMOList,
                             batch = NULL,
                             col1 = "#440154",
                             col2 = "#fde725") {
  # Check if smallRNA data or annotation is available
  if (is.null(getRawData(objMOList, "smallRNAseq"))) {
    stop("No small RNA data available.")
  } else if (is.null(objMOList$annoSncRNA)) {
    stop("No annotation for small RNA data. Please run annotateSmallRNA first.")
  } else {
    # Continue
  }

  # Retrieve the raw expression matrix for small RNA annotation
  dfSncRNA <- objMOList %>%
    getRawData("smallRNAseq") %>%
    as.data.frame()
  sncAnno <- objMOList$annoSncRNA
  if (all(sncAnno == HUMAN)) {
    dfSncRNA <- annoSncList(dfSncRNA,
      annoList = SNCANNOLIST_HSAPIENS
    )
  } else if (is.list(sncAnno)) {
    dfSncRNA <- annoSncList(dfSncRNA, sncAnno)
  } else {
    stop("Incompatible annotation for small RNA. Please check the annotation.")
  }

  # Separate the expression matrix for each type of small RNA
  sncTypes <- unique(dfSncRNA$type)
  if (NOANNO %in% sncTypes) {
    numUnannotated <- sum(dfSncRNA$type == NOANNO)
    warning(paste0("Unannotated small RNAs: ", numUnannotated))
    warning("Some small RNAs are not annotated. Unannotated small RNAs are
    excluded from the analysis.")
    dfSncRNA <- dfSncRNA %>% dplyr::filter(type != NOANNO)
    sncTypes <- sncTypes[sncTypes != NOANNO]

    # Check if there is any small RNA type left
    if (length(sncTypes) == 0) {
      # should not reach here if the users do not temper the annotation
      stop("Error validating small RNA types. Please check the annotation.")
    } else {
      # Continue
    }
  } else {
    # Continue
  }

  # Generate a list of raw expression matrices
  matrixSncRNAList <- list()
  # Separate by type and remove the type variable for a numeric matrix
  for (sncType in sncTypes) {
    matrixSncRNAList[[sncType]] <- dfSncRNA %>%
      dplyr::filter(type == sncType) %>%
      dplyr::select(-type) %>%
      as.matrix()
  }

  # Perform PCA calculation and generate the plots
  pcaPlotList <- list()
  for (sncType in sncTypes) {
    cat("Performing PCA analysis for ", sncType, "...\n", sep = "")
    pcaPlotList[[sncType]] <- countPCA(
      matrix = matrixSncRNAList[[sncType]],
      groupBy = objMOList@smallRNAseqSamples$groupBy,
      title = sncType
    )
  }
  cat("Finished estimating PCA for existing small RNA types.\n")
  cat("Individual PCA plots can be assessed by pcaPlotList$<small RNA type>.\n")

  rm(dfSncRNA, matrixSncRNAList) # free up memory from big objects

  return(pcaPlotList)
}


#' Peak coverage plot for differential accessible regions
#'
#' @aliases plotATACCoverage
#'
#' @description This function generates a peak coverage plot for differential
#'              accessible regions based on the ATACseq data. This is a wrapper
#'              function for the covplot function in the ChIPseeker package.
#'
#' @param objMOList An MOList object containing the differential accessible
#'                  regions
#' @param title The title for the plot, default is "ATAC Peaks over Chromosomes"
#'
#' @return A ggplot object
#'
#' @importFrom ChIPseeker covplot
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#'
#' @export
#'
#' @references
#' \insertRef{yu2015chipseeker}{IntegraTRN}
#'
#' \insertRef{lawrence2013software}{IntegraTRN}
#'
#' @examples
#' \dontrun{
#' # Use the package-provided example data
#' data(expMOList)
#'
#' # Annotates the peaks in the example data and performs motif enrichment
#' # Using human annotations as an example
#' library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' library("org.Hs.eg.db")
#' library("BSgenome.Hsapiens.UCSC.hg38")
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' annoDb <- "org.Hs.eg.db"
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38
#'
#' # Load PWMs from JASPAR
#' data("jasparVertebratePWM")
#'
#' # Performs the annotation and motif enrichment
#' expMOList <- annotateATACPeaksMotif(expMOList,
#'   tssRegion = c(-3000, 3000),
#'   TxDb = txdb,
#'   annoDb = annoDb,
#'   bsgenome = bsgenome,
#'   pwmL = jasparVertebratePWM
#' )
#'
#' # Plotting the coverage plot
#' plotATACCoverage(expMOList)
#' }
#'
plotATACCoverage <- function(objMOList, title = "ATAC Peaks over Chromosomes") {
  if (is.null(objMOList$DEATAC)) {
    stop("No differential accessible regions. Please run diffOmics() first.
  See ?diffOmics for details.")
  } else {
    # Continue
  }
  if (inherits(objMOList$DEATAC, "PEAKTag")) {
    # A PEAKTag object
    peakGR <- asGRanges(objMOList$DEATAC)
  } else {
    # A DETag object
    peakGR <- GenomicRanges::makeGRangesFromDataFrame(exportDE(objMOList$DEATAC))
  }
  covPlot <- ChIPseeker::covplot(peakGR, title = title)
  return(covPlot)
}


#' Plotting the annotation Pie chart of differential accessible regions
#'
#' @aliases plotATACAnno
#'
#' @description This function generates a pie chart for the annotation of
#'              differential accessible regions based on the ATACseq data.
#'              This is a wrapper function for the plotAnnoPie function in the
#'              ChIPseeker package.
#'
#' @param objMOList An MOList object containing the differential accessible
#'                  regions
#'
#' @return A ggplot object
#'
#' @importFrom ChIPseeker plotAnnoPie
#'
#' @export
#'
#' @references
#' \insertRef{yu2015chipseeker}{IntegraTRN}
#'
#' @examples
#' \dontrun{
#' # Use the package-provided example data
#' data(expMOList)
#'
#' # Annotates the peaks in the example data and performs motif enrichment
#' # Using human annotations as an example
#' library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' library("org.Hs.eg.db")
#' library("BSgenome.Hsapiens.UCSC.hg38")
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' annoDb <- "org.Hs.eg.db"
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38
#'
#' # Load PWMs from JASPAR
#' data("jasparVertebratePWM")
#'
#' # Performs the annotation and motif enrichment
#' expMOList <- annotateATACPeaksMotif(expMOList,
#'   tssRegion = c(-3000, 3000),
#'   TxDb = txdb,
#'   annoDb = annoDb,
#'   bsgenome = bsgenome,
#'   pwmL = jasparVertebratePWM
#' )
#'
#' # Plotting the annotation pie chart
#' plotATACAnno(expMOList)
#' }
#'
plotATACAnno <- function(objMOList) {
  if (is.null(objMOList$DEATAC) ||
    !inherits(objMOList$DEATAC, "PEAKTag") ||
    is.null(objMOList$DEATAC@annotatedPeaks)) {
    stop("No peak annotation found.")
  } else {
    # Continue
  }
  csAnno <- objMOList$DEATAC@annotatedPeaks
  ChIPseeker::plotAnnoPie(csAnno)
}


#' Motif heatmap of differential accessible regions of enriched motifs
#'
#' @aliases plotATACMotifHeatmap
#'
#' @description This function generates a combined profile and heatmap of
#'              differential accessible regions of enriched motifs based on the
#'              ATACseq data. This is a wrapper function for the
#'              plotMotifHeatmaps function in the ChIPseeker package.
#'
#' @details The motif enrichment analysis can be exceptionally stringent for
#'          ATACseq data. The users can specify the cutoff for the adjusted
#'          p-value or the regular p-value. The users can also specify the
#'          log2 enrichment fold change cutoff. Please use this function to
#'          explore the enriched motifs in the differential accessible regions
#'          and decide the best combination of cutoffs for constructing the
#'          transcriptional regulatory network. If pValue is specified,
#'          pValueAdj will be ignored.
#'
#' @param objMOList An MOList object containing the differential accessible
#'                  regions
#' @param pValueAdj The cutoff for adjusted p-value, default is 0.05
#' @param pValue The cutoff for p-value. If pValue is specified, pValueAdj
#'               will be ignored
#' @param log2FEnrich The cutoff for log2 fold enrichment, default is NULL
#'
#' @return A ComplexHeatmap object
#'
#' @importFrom monaLisa plotMotifHeatmaps
#' @importFrom SummarizedExperiment assay
#'
#' @export
#'
#' @references
#' \insertRef{machlab2022monalisa}{IntegraTRN}
#'
#' \insertRef{SummarizedExperiment}{IntegraTRN}
#'
#' \insertRef{gu2016complex}{IntegraTRN}
#'
#'
#' @examples
#' \dontrun{
#' # Use the package-provided example data
#' data(expMOList)
#'
#' # Annotates the peaks in the example data and performs motif enrichment
#' # Using human annotations as an example
#' library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' library("org.Hs.eg.db")
#' library("BSgenome.Hsapiens.UCSC.hg38")
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' annoDb <- "org.Hs.eg.db"
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38
#'
#' # Load PWMs from JASPAR
#' data("jasparVertebratePWM")
#'
#' # Performs the annotation and motif enrichment
#' expMOList <- annotateATACPeaksMotif(expMOList,
#'   tssRegion = c(-3000, 3000),
#'   TxDb = txdb,
#'   annoDb = annoDb,
#'   bsgenome = bsgenome,
#'   pwmL = jasparVertebratePWM
#' )
#'
#' # Example 1: Plot the motif heatmap by default parameters
#' plotATACMotifHeatmap(expMOList)
#'
#' # Example 2: Plot the motif heatmap with unadjusted p-value cutoff of 0.01
#' plotATACMotifHeatmap(expMOList, pValue = 0.01)
#'
#' # Example 3: Plot the motif heatmap with log2 fold enrichment cutoff of 1
#' plotATACMotifHeatmap(expMOList, pValue = 0.01, log2FEnrich = 1)
#' }
#'
plotATACMotifHeatmap <- function(objMOList,
                                 pValueAdj = 0.05,
                                 pValue = NULL,
                                 log2FEnrich = NULL) {
  if (is.null(objMOList$DEATAC) ||
    !inherits(objMOList$DEATAC, "PEAKTag") ||
    is.null(objMOList$DEATAC$motifEnrichment)) {
    stop("No motif enrichment analysis found.")
  } else {
    # Do nothing
  }

  # Select enriched motifs
  enrichedMotifs <- objMOList$DEATAC$motifEnrichment
  sel <- selectedMotifs(
    enrichedMotifs = enrichedMotifs,
    pValueAdj = pValueAdj,
    pValue = pValue,
    log2FEnrich = log2FEnrich
  )

  # Define values to plot
  if (is.null(pValue)) {
    pramToPlot <- c("log2enr", "negLog10Padj")
    maxSig <- max((enrichedMotifs[sel, ])@assays@data@listData$negLog10Padj)
  } else {
    pramToPlot <- c("log2enr", "negLog10P")
    maxSig <- max((enrichedMotifs[sel, ])@assays@data@listData$negLog10P)
  }

  # Styling parameter
  maxEnr <- max((enrichedMotifs[sel, ])@assays@data@listData$log2enr)

  # Plot the heatmap
  monaLisa::plotMotifHeatmaps(
    x = enrichedMotifs[sel, ],
    which.plots = pramToPlot,
    width = 1.8,
    cluster = TRUE,
    maxEnr = maxEnr,
    maxSig = maxSig,
    show_seqlogo = TRUE
  )
}


# [END]
