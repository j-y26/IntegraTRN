# Purpose: Plotting differentially expressed omics data
# Author: Jielin Yang
# Date: 2023-11-02
# Version: 1.0
# Bugs and Issues: None


#' Annotate by expression
#'
#' @description This function annotates whether a gene is up-regulated,
#'              down-regulated, or not differentially expressed
#'
#' @param deg A data frame containing the differential expression analysis
#'            results, must follow the described format in DETag-class
#' @param log2FC The cutoff for log2 fold change
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
#' @importFrom ggplot2 ylab theme element_blank
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
    ggplot2::geom_vline(xintercept = c(-log2FC, log2FC), linetype = "dashed")

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
#' @importFrom ggplot2 ylab theme element_blank
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
annoSncDEList <- function(deg, annoList) {
  deg <- deg %>% dplyr::mutate(type = dplyr::case_when(
    rownames(deg) %in% annoList$miRNA ~ "miRNA",
    rownames(deg) %in% annoList$piRNA ~ "piRNA",
    rownames(deg) %in% annoList$snoRNA ~ "snoRNA",
    rownames(deg) %in% annoList$snRNA ~ "snRNA",
    rownames(deg) %in% annoList$tRNA ~ "tRNA",
    rownames(deg) %in% annoList$circRNA ~ "circRNA",
    TRUE ~ "Unannotated"
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
    degSmallRNAseq <- annoSncDEList(degSmallRNAseq,
      annoList = list(
        miRNA = miRNAHsapiens,
        piRNA = piRNAHsapiens,
        snoRNA = snoRNAHsapiens,
        snRNA = snRNAHsapiens,
        tRNA = tRNAHsapiens,
        circRNA = circRNAHsapiens
      )
    )
  } else if (is.list(sncAnno)) {
    degSmallRNAseq <- annoSncDEList(degSmallRNAseq, sncAnno)
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
