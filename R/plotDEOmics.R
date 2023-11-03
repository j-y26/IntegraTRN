# Purpose: Plotting differentially expressed omics data
# Author: Jielin Yang
# Date: 2023-11-02
# Version: 1.0
# Bugs and Issues: None


# Define the global variables
NOANNO <- "Unannotated"


#' Annotate by expression
#'
#' @description This function annotates whether a gene is up-regulated,
#'              down-regulated, or not differentially expressed
#'
#' @param deg A data frame containing the differential expression analysis
#'            results, must follow the described format in DETag-class
#' @param log2FC The cutoff for log2 fold change, a positive number
#' @param adjP The cutoff for adjusted p-value
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
#' @importFrom ggrepel geom_label_repel
#'
plotVolcano <- function(deg,
                        log2FC = 1,
                        adjP = 0.05,
                        geneList = NULL,
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

  # Label the genes in the geneList
  labelData <- deg[rownames(deg) %in% geneList, ]
  if (!is.null(geneList)) {
    vPlot <- vPlot + ggrepel::geom_label_repel(
      data = labelData,
      mapping = ggplot2::aes(logFC, -log(padj, 10),
        label = rownames(labelData)
      ),
      size = 3,
      color = "black",
      nudge_x = 0.1,
      nudge_y = 0.1
    )
  } else {
    # Do nothing
  }

  # Add the title
  if (!is.null(title)) {
    vPlot <- vPlot + ggplot2::ggtitle(title)
  }

  return(vPlot)
}


#' Volcano plot for visualizing differentially expressed mRNA
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
#' @param log2FC The cutoff for log2 fold change. Default is 1.
#' @param adjP The cutoff for adjusted p-value. Default is 0.05.
#' @param geneList A vector of genes to highlight in the plot. Default is NULL.
#' @param upColor The color for up-regulated genes. Default is "firebrick3".
#' @param downColor The color for down-regulated genes. Default is "dodgerblue3"
#' @param title The title for the plot. Default is NULL.
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 aes scale_color_manual guides guide_legend theme_bw xlab
#' @importFrom ggplot2 ylab theme element_blank element_text element_line
#' @importFrom ggplot2 ggplot geom_point ggtitle geom_text
#' @importFrom dplyr mutate case_when
#' @importFrom ggrepel geom_label_repel
#'
#' @export
#'
#' @examples
#' # Suppose that we have an MOList object called objMOList, which contains the
#' # differential expression results for mRNA as an element named DERNAseq.
#'
#' # Example 1: Generate the volcano plot by default parameters
#' \dontrun{
#' plotVolcanoRNA(objMOList)
#' }
#'
#' # Example 2: Generate the volcano plot with custom parameters
#' \dontrun{
#' plotVolcanoRNA(objMOList,
#'   log2FC = 2,
#'   adjP = 0.01,
#'   geneList = c("MYH7B", "NELFCD", "KDM8"),
#'   upColor = "purple",
#'   downColor = "green",
#'   title = "Volcano plot for DE mRNA"
#' )
#' }
#'
plotVolcanoRNA <- function(objMOList,
                           log2FC = 1,
                           adjP = 0.05,
                           geneList = NULL,
                           upColor = "firebrick3",
                           downColor = "dodgerblue3",
                           title = NULL) {
  if (is.null(objMOList$DERNAseq)) {
    stop("No differential expression results for mRNA. Please run diffOmics()
  first. See ?diffOmics for details.")
  } else {
    # Continue
  }

  # Retrieve the differential expression results for mRNA
  degRNAseq <- objMOList$DERNAseq %>% exportDE()

  # Annotate expression
  degRNAseq <- annoExpr(degRNAseq, log2FC, adjP)

  # Generate the base volcano plot
  vPlot <- plotVolcano(degRNAseq, log2FC, adjP, geneList, title)

  # Color the up- and down-regulated genes
  vPlot <- vPlot +
    ggplot2::geom_point(ggplot2::aes(color = expr), size = 4 / 5) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes =
        list(size = 2.5)
    )) +
    ggplot2::scale_color_manual(values = c(
      "Up-regulated" = upColor,
      "Down-regulated" = downColor,
      "Not DE" = "grey50"
    ))

  # Annotate the numbers of up- and down-regulated genes
  upNum <- sum(degRNAseq$logFC >= log2FC & degRNAseq$padj < adjP)
  downNum <- sum(degRNAseq$logFC <= -log2FC & degRNAseq$padj < adjP)
  vPlot <- vPlot + ggplot2::geom_text(
    ggplot2::aes(x = max(degRNAseq$logFC), y = -max(log(degRNAseq$padj, 10))),
    label = upNum,
    color = upColor,
    hjust = 1,
    vjust = 1
  ) + ggplot2::geom_text(
    ggplot2::aes(x = min(degRNAseq$logFC), y = -max(log(degRNAseq$padj, 10))),
    label = downNum,
    color = downColor,
    hjust = 0,
    vjust = 1
  )

  return(vPlot)
}


#' Annotate the type of small RNA in the differential expression results
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
#' @param log2FC The cutoff for log2 fold change. Default is 1.
#' @param adjP The cutoff for adjusted p-value. Default is 0.05.
#' @param geneList A vector of genes to highlight in the plot.
#' @param colScheme A RColorBrewer color scheme for color-coding each type of
#'                  small RNAs. Default is "BuPu".
#' @param title The title for the plot. Default is NULL.
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 aes scale_color_manual guides guide_legend theme_bw xlab
#' @importFrom ggplot2 ylab theme element_blank
#' @importFrom ggplot2 ggplot geom_point ggtitle geom_text
#' @importFrom dplyr mutate case_when
#' @importFrom ggrepel geom_label_repel
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
#'
#' @examples
#' # Suppose that we have an MOList object called objMOList, which contains the
#' # differential expression results for small RNAs as an element named
#' DESmallRNAseq
#'
#' # Example 1: Generate the volcano plot by default parameters
#' \dontrun{
#' plotVolcanoSmallRNA(objMOList)
#' }
#'
#' # Example 2: Generate the volcano plot with custom parameters
#' \dontrun{
#' plotVolcanoSmallRNA(objMOList,
#'   log2FC = 2,
#'   adjP = 0.01,
#'   geneList = c("hsa-miR-1-3p", "hsa-miR-2-3p"),
#'   colScheme = "magma",
#'   title = "Volcano plot for DE small RNA"
#' )
#' }
#'
plotVolcanoSmallRNA <- function(objMOList,
                                log2FC = 1,
                                adjP = 0.05,
                                geneList = NULL,
                                colScheme = "BuPu",
                                title = NULL) {
  if (is.null(objMOList$DESmallRNAseq)) {
    stop("No differential expression results for small RNA. Please run
  diffOmics() first. See ?diffOmics for details.")
  } else if (is.null(objMOList$annoSncRNA)) {
    stop("No annotations for small RNA. Please run annotateSmallRNA() first.")
  } else {
    # Continue
  }

  # Retrieve the differential expression results for small RNA
  degSmallRNAseq <- objMOList$DESmallRNAseq %>% exportDE()

  # Annotate the type of small RNA
  sncAnno <- objMOList$annoSncRNA
  if (sncAnno == HUMAN) {
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

  # Generate the base volcano plot
  vPlot <- plotVolcano(degSmallRNAseq, log2FC, adjP, geneList, title)


  # Generate a set of colors for each type of small RNA
  sncColors <- RColorBrewer::brewer.pal(length(sncTypes), colScheme)
  sncColorList <- setNames(sncColors, sncTypes)


  # Color each type of small RNA
  vPlot <- vPlot +
    ggplot2::geom_point(ggplot2::aes(color = type), size = 4 / 5) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes =
        list(size = 2.5)
    )) +
    ggplot2::scale_color_manual(values = c(
      sncColorList,
      "Not DE" = "grey50"
    ))
  return(vPlot)
}


#' Principle Component Analysis (PCA) for count-based omics data
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
#' @importFrom ggrepel geom_label_repel
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq
#' @importFrom DESeq2 varianceStabilizingTransformation plotPCA
#'
#' @return A ggplot object
#'
#' @export
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
    ggplot2::coord_fixed() +
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
  if (length(unique(groupBy)) > 2) {
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
#' @importFrom ggrepel geom_label_repel
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq
#' @importFrom DESeq2 varianceStabilizingTransformation plotPCA
#'
#' @export
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
  if (sncAnno == HUMAN) {
    dfSncRNA <- annoSncList(dfSncRNA,
      annoList = SNCANNOLIST_HSAPIENS
    )
  } else if (is.list(sncAnno)) {
    dfSncRNA <- annoSncList(degSmallRNAseq, sncAnno)
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
      title = sncType
    )
  }
  cat("Finished estimating PCA for existing small RNA types.\n")
  cat("Individual PCA plots can be assessed by pcaPlotList$<small RNA type>.\n")

  rm(dfSncRNA, matrixSncRNAList) # free up memory from big objects

  return(pcaPlotList)
}


# [END]
