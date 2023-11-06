# Purpose: Processing peak files and motif enrichment analysis
# Author: Jielin Yang
# Date: 2023-11-05
# Version: 1.0
# Bugs and Issues: None


#' Merge peaks with overlapping genomic coordinates and filter out highly
#' overlapped peaks between conditions
#'
#' @keywords internal
#'
#' @param peaks A data frame containing the peak information, in BED format
#'
#' @return A GRanges object containing the merged peaks
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame reduce
#'
mergePeaks <- function(peaks) {
  # Generate a GRanges object from the peak data frame
  # The GRangesFromDataFrame function requires "start" and "end" as column names
  # No strandedness in ATACseq
  peakGR <- GenomicRanges::makeGRangesFromDataFrame(peaks,
    ignore.strand = TRUE
  )

  # Merge any overlapping peaks into a single peak
  # Reduce internally uses the inter-range transformation from the IRanges
  # package
  peakGR <- GenomicRanges::reduce(peakGR)
  return(peakGR)
}


#' Processing user input peak files
#'
#' @keywords internal
#'
#' @description This function processes user input peak files to check for
#'              overlaps in genomic coordinates. Although this process is not
#'              necessary for high-quality peak files, it improves the
#'              robustness of the subsequent motif enrichment analysis and
#'              peak annotation. Any peaks that are found to be overlapping
#'              with a minimum of 1 bp overlap will be merged into a single
#'              peak, for peaks from a single condition. Since the two input
#'              peak files that already represents differential accessible
#'              regions, we do not merge peaks across conditions. Instead,
#'              peaks with overlapping genomic coordinates that covers more than
#'              50% of the peak region will be removed.
#'
#' @param objMOList An object of class MOList
#'
#' @return An object of class MOList, with a DETag added containing a master
#'         peak set that annotates where the peaks come from
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame reduce findOverlaps
#' @importFrom IRanges width
#' @importFrom Repitools annoGR2DF
#' @importFrom S4Vectors from to
#'
processPeakOverlap <- function(objMOList) {
  # Check if the ATACseq peak files exist
  peaksCond1 <- getRawData(objMOList, "ATAC")[[1]]
  peaksCond2 <- getRawData(objMOList, "ATAC")[[2]]
  if (is.null(peaksCond1) || is.null(peaksCond2)) {
    return(objMOList)
  }

  # Respectively check for peak overlaps within each condition
  peakGR1 <- mergePeaks(peaksCond1)
  peakGR2 <- mergePeaks(peaksCond2)

  # Remove overlapping peaks across conditions, with having at least 50% overlap
  # in any peak considered as overlapping
  overlaps <- GenomicRanges::findOverlaps(peakGR1, peakGR2)
  fromP1 <- S4Vectors::from(overlaps) # index of overlapped peaks in peakGR1
  fromP2 <- S4Vectors::to(overlaps) # index of pverlapped peaks in peakGR2
  overlappedWidth <- IRanges::width(GenomicRanges::pintersect(
    peakGR1[fromP1],
    peakGR2[fromP2]
  ))
  # Calculate the minimum peak width for peaks that overlap
  width1 <- IRanges::width(peakGR1[fromP1])
  width2 <- IRanges::width(peakGR2[fromP2])

  peakGR1 <- peakGR1[-fromP1[overlappedWidth / width1 > 0.5]]
  peakGR2 <- peakGR2[-fromP2[overlappedWidth / width2 > 0.5]]

  # Convert the GRanges objects back to data frames
  peakDF1 <- Repitools::annoGR2DF(peakGR1)
  peakDF2 <- Repitools::annoGR2DF(peakGR2)

  # Make a master peak set
  peakDF1 <- peakDF1 %>%
    dplyr::select(-width) %>%
    dplyr::mutate(Condition = "-")
  peakDF2 <- peakDF2 %>%
    dplyr::select(-width) %>%
    dplyr::mutate(Condition = "+")
  masterPeaks <- dplyr::bind_rows(peakDF1, peakDF2)

  atacDETag <- DETag(masterPeaks, ATAC_GRANGE)

  # Update the MOList object
  objMOList$DEATAC <- atacDETag
  return(objMOList)
}


#' Annotate peaks with genomic features
#'
#' @description This function annotates the ATACseq peaks with genomic features
#'              using the ChIPseeker package.
#'
#' @details The annotation is performed using the TxDb object and annotation
#'          database specified by the user. Depending on the style of genomic
#'          coordinate representation (with or without the "chr" prefix for the
#'          chromosome name, or the use of "chrM", "M", "MT" for mitochondrial
#'          DNA), the TxDb object may need to be adjusted to match the style
#'          used in the peak file. See the vignette of the ChIPseeker package
#'          for more details. By default, we recommend using the UCSC style
#'          genomic coordinates, and the use of the
#'          TxDb.Hsapiens.UCSC.hg38.knownGene and org.Hs.eg.db packages for
#'          human samples, and the TxDb.Mmusculus.UCSC.mm10.knownGene and
#'          org.Mm.eg.db packages for mouse samples.
#'
#' @param objMOList An object of class MOList
#' @param tssRegion The region around the TSS to annotate with
#'                  (default: +/- 3000)
#' @param TxDB The TxDb object to use for annotation, using one of the TxDb
#'             packages
#' @param annoDb The annotation database to use for annotation, using one of
#'               the org.*.eg.db packages. Must be a valid string for the name
#'               of the package. For details, see the vignette of the ChIPseeker
#'               package for custom annotation databases.
#'
#' @return A csAnno object containing the annotated peaks, as defined by the
#'         ChIPseeker package. In particular, the original "Condition"
#'         annotation is preserved.
#'
#' @importFrom ChIPseeker annotatePeak
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#'
annotatePeaks <- function(objMOList,
                          tssRegion = c(-3000, 3000),
                          TxDb,
                          annoDb = "org.Hs.eg.db") {
  # Retrieve the ATACseq peaks
  peaks <- exportDE(objMOList$DEATAC)

  # Generate a GRanges object from the peak data frame
  peakGR <- GenomicRanges::makeGRangesFromDataFrame(peaks,
    keep.extra.columns = TRUE
  )

  # Annotate the peaks with genomic features
  # The annotation is performed using the TxDb object and annotation database
  # specified by the user
  peakAnno <- ChIPseeker::annotatePeak(
    peakGR,
    tssRegion = tssRegion,
    TxDb = TxDb,
    annoDb = annoDb
  )
  # Generate a PEAKTag object that replaces original DETag object
  peakTag <- PEAKTag(objMOList$DEATAC, peakAnno, TxDb, annoDb)
  objMOList$DEATAC <- peakTag
  return(objMOList)
}
