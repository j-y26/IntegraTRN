# Purpose: Plotting differentially expressed omics data
# Author: Jielin Yang
# Date: 2023-11-02
# Version: 1.0
# Bugs and Issues: None



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
#' @param upColor The color for up-regulated genes.
#' @param downColor The color for down-regulated genes.
#' @param title The title for the plot.
#' 
#' @return A ggplot object
#' 
#' @importFrom ggplot2 aes color scale_color_manual guides theme_bw xlab ylab
#' @importFrom ggplot2 ggplot geom_point ggtitle theme element_blank
#' @importFrom dplyr mutate case_when
#' @importFrom stats log
#' 
#' 
plotVolcano <- function(deg, log2FC, adjP, upColor, downColor, title) {
  # Annotate the data for plotting
  deg <- deg %>%
    dplyr::mutate(expr = dplyr::case_when(
      logFC >= log2FC & padj < adjP ~ "Up-regulated",
      logFC <= -log2FC & padj < adjP ~ "Down-regulated",
      TRUE ~ "Not DE"))
  
  # Generate the plot
  vPlot <- ggplot2::ggplot(deg, ggplot2::aes(x = logFC, y = -log(padj, 10))) +
    ggplot2::geom_point(ggplot2::aes(color = expr), size = 3/5) +
    ggplot2::scale_color_manual(values = c("Up-regulated" = upColor, 
                                           "Down-regulated" = downColor,
                                           "Not DE" = "grey50")) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes =
                                                      list(size = 2.5))) +
    ggplot2::theme_bw() +
    ggplot2::xlab(expression("log"[2]*"FC")) + 
    ggplot2::ylab(expression("-log"[10]*"adj.p-value")) + 
    ggplot2::theme(legend.title = ggplot2::element_blank())
  
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
#' @importFrom ggplot2 aes color scale_color_manual guides theme_bw xlab ylab
#' @importFrom ggplot2 ggplot geom_point ggtitle geom_text theme element_blank
#' @importFrom dplyr mutate case_when
#' @importFrom stats log
#'
#' @export
#'
#' @examples
#' 
#'
#'
plotVolcanoRNA <- function(objMOList,
                           log2FC = 1,
                           adjP = 0.05,
                           geneList = NULL,
                           upColor = "red",
                           downColor = "blue",
                           title = NULL) {

  if (is.null(objMOList$DERNAseq)) {
    stop("No differential expression results for mRNA. Please run diffOmics()
         first. See ?diffOmics for details.")
  }

  # Retrieve the differential expression results for mRNA
  degRNAseq <- objMOList$DERNAseq %>% exportDE()

  # Generate the base volcano plot
  vPlot <- plotVolcano(degRNAseq, log2FC, adjP, upColor, downColor, title)

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
