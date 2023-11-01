# Purpose: Performing differential accessibility and motif analysis on ATAC-seq 
#          peaks
# Author: Jielin Yang
# Date: 2023-10-31
# Version: 1.0
# Bugs and Issues: None


#' Identify differential peaks from two ATAC-seq peak sets in BED format
#'
#' @description This function identifies differential peaks from two ATAC-seq
#'              peak sets in BED format. The function generates a global peak
#'              set that accounts for regions that are present in one but not
#'              both BED files. The function then performs annotates the peaks
#'              with which condition they are present. A peak that is annotated
#'              as "+" is identified to present only in condition 2, and vice
#'              versa for "-".
#'
#' @param objMOList an object of class MOList
#'
#' @return A data frame containing the differential peaks, in BED4 format
#' \describe{
#' \item{chrom}{Chromosome name}
#' \item{chromStart}{Start position of the peak}
#' \item{chromEnd}{End position of the peak}
#' \item{name}{Annotation of the peak about its condition, either "+" or "-"}
#'
#'
diffPeaks <- function(objMOList) {
  # Determine whether ATACseq data exists to decide downstream analysis
  peaksCond1 <- getRawData(objMOList, "ATAC")[["peaksCond1"]]
  peaksCond2 <- getRawData(objMOList, "ATAC")[["peaksCond2"]]
  if (is.null(peaksCond1) || is.null(peaksCond2)) {
    return(objMOList)
  }

  

  # Identify overlapping peaks
  peakGRanges <- GenomicRanges::reduce(
    GenomicRanges::union(peakGRanges1, peakGRanges2)
  )
}
