# Purpose: Processing peak files and motif enrichment analysis
# Author: Jielin Yang
# Date: 2023-11-01
# Version: 1.0
# Bugs and Issues: None


#' Merge peaks with overlapping genomic coordinates
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
  colnames(peaks)[2:3] <- c("start", "end")
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
#'              regions, we do not merge peaks across conditions.
#'
#' @param objMOList An object of class MOList
#'
#' @return An object of class MOList, with a DETag added containing a master
#'         peak set that annotates where the peaks come from
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame reduce
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

  # Master peak set
  masterPeaks <- data.frame(
    seqnames = peakGR1@seqnames,
    start = peakGR1@ranges@start,
    end = start + peakGR1@ranges@width - 1,
    Cond = "Cond1"
  )
  masterPeaks <- rbind(
    masterPeaks,
    data.frame(
      seqnames = peakGR2@seqnames,
      start = peakGR2@ranges@start,
      end = start + peakGR2@ranges@width - 1,
      Cond = "Cond2"
    )
  )
  atacDETag <- DETag(masterPeaks, ATAC_GRANGE)

  # Update the MOList object
  objMOList$DEATAC <- atacDETag
  return(objMOList)
}
