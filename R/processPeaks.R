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
#' @references
#' \insertRef{lawrence2013software}{IntegraTRN}
#'
mergePeaks <- function(peaks) {
  # Generate a GRanges object from the peak data frame
  # The GRangesFromDataFrame function requires "start" and "end" as column names
  # No strandedness in ATACseq
  peakGR <- GenomicRanges::makeGRangesFromDataFrame(peaks,
    ignore.strand = TRUE
  )

  # Merge any overlapping peaks into a single peak
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
#' @importFrom GenomicRanges makeGRangesFromDataFrame reduce findOverlaps
#'
#' @references
#' \insertRef{lawrence2013software}{IntegraTRN}
#'
processPeakOverlap <- function(objMOList) {
  # Check if the ATACseq peak files exist
  peaksCond1 <- getRawData(objMOList, "ATAC")[[1]]
  peaksCond2 <- getRawData(objMOList, "ATAC")[[2]]
  if (is.null(peaksCond1) || is.null(peaksCond2)) {
    return(objMOList)
  }

  cat("Processing peaks...\n")

  # Respectively check for peak overlaps within each condition
  peakGR1 <- mergePeaks(peaksCond1)
  peakGR2 <- mergePeaks(peaksCond2)

  # Convert the GRanges objects back to data frames
  peakDF1 <- as.data.frame(peakGR1)[, 1:4]
  colnames(peakDF1)[1:3] <- CHROMINFO
  peakDF2 <- as.data.frame(peakGR2)[, 1:4]
  colnames(peakDF2)[1:3] <- CHROMINFO

  # Make a master peak set
  peakDF1 <- peakDF1 %>%
    dplyr::select(-width) %>%
    dplyr::mutate(Condition = "-")
  peakDF2 <- peakDF2 %>%
    dplyr::select(-width) %>%
    dplyr::mutate(Condition = "+")
  masterPeaks <- dplyr::bind_rows(peakDF1, peakDF2)

  atacDETag <- DETag(masterPeaks, ATAC_GRANGE)

  # remove unnecessary objects
  rm(peakGR1, peakGR2, peakDF1, peakDF2, masterPeaks)

  # Update the MOList object
  objMOList$DEATAC <- atacDETag
  return(objMOList)
}


#' Annotate peaks with genomic features
#'
#' @keywords internal
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
#'          for more details. By default, it is recommended using the UCSC style
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
#' @return A MOList object containing a csAnno object of annotated peaks, as
#'         defined by the
#'         ChIPseeker package. In particular, the original "Condition"
#'         annotation is preserved.
#'
#' @importFrom ChIPseeker annotatePeak
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#'
#' @references
#' \insertRef{yu2015chipseeker}{IntegraTRN}
#'
#' \insertRef{lawrence2013software}{IntegraTRN}
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


#' Perform motif enrichment analysis on peaks with binary conditions
#'
#' @keywords internal
#'
#' @description This function performs motif enrichment analysis on peaks with
#'              binary conditions. The motif enrichment analysis is performed
#'              using the monaLisa package.
#'
#' @param objMOList An object of class MOList
#' @param bsgenome The BSgenome object that contains the genome sequence
#'                 Read the Bioconductor
#'                 vignette for the BSgenome package for more details.
#' @param pwmL A PWMatrixList object containing the position-weight matrices
#'             (PWMs) to use for motif enrichment analysis. The PWMs can be
#'             generated using the getMatrixSet() function from the TFBSTools
#'             package. See the vignette of the TFBSTools package for more
#'             details.
#' @param fixedWidth The fixed width to adjust the peak size for motif if
#'                   necessary (default: 500). The default width follows an
#'                   optimal width peak selection strategy defined by
#'                   \insertCite{grandi2022chromatin;textual}{IntegraTRN}.
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame resize width
#' @importFrom monaLisa calcBinnedMotifEnrR
#' @importFrom Biostrings getSeq
#'
#' @references
#' \insertRef{BSgenome}{IntegraTRN}
#'
#' \insertRef{tan2016tfbstools}{IntegraTRN}
#'
#' \insertRef{grandi2022chromatin}{IntegraTRN}
#'
#' \insertRef{machlab2022monalisa}{IntegraTRN}
#'
#' \insertRef{Biostrings}{IntegraTRN}
#'
enrichMotifs <- function(objMOList, bsgenome, pwmL, fixedWidth = 500) {
  if (is.null(objMOList$DEATAC)) {
    # No ATACseq peaks
    stop("No differential accessible regions. Please run diffOmics() first.
    See ?diffOmics for details.")
  } else if (!inherits(objMOList$DEATAC, "PEAKTag")) {
    # The ATACseq peaks are not annotated
    stop("The ATACseq peaks are not annotated. Please run annotatePeaks()
    first. See ?annotatePeaks for details.")
  } else {
    # Retrieve the ATACseq peaks
    peaks <- exportDE(objMOList$DEATAC)
    cond <- peaks$Condition
    bins <- factor(cond)
  }

  # Check the validity of the PWMs
  if (!inherits(pwmL, "PWMatrixList")) {
    stop("The PWMs must be a PWMatrixList object. See ?enrichMotifs for
    details.")
  } else {
    # Continue
  }

  cat("Performing motif enrichment analysis...this may take a while...\n")

  # Generate GRanges mobject
  peakGR <- GenomicRanges::makeGRangesFromDataFrame(peaks,
    keep.extra.columns = TRUE
  )

  # First check if the peak size is the same
  peakWidth <- GenomicRanges::width(peakGR)
  if (length(unique(peakWidth)) > 1) {
    # The peak size is not the same, need to fix the peak size for motif
    # enrichment analysis
    peakGR <- GenomicRanges::resize(peakGR, width = fixedWidth, fix = "center")
  } else {
    # The peak size is the same
    # No need to fix the peak size
  }

  # Generate a DNAStringSet object for motif enrichment
  peakSeq <- Biostrings::getSeq(bsgenome, peakGR)

  # Perform motif enrichment analysis on binary conditions
  enrichedMotifs <- monaLisa::calcBinnedMotifEnrR(
    seqs = peakSeq,
    bins = bins,
    pwmL = pwmL,
  )

  cat("Done.\n")

  # Remove large objects
  rm(peakGR, peakSeq)

  # Append the motif enrichment results to the MOList object
  objMOList$DEATAC$motifEnrichment <- enrichedMotifs
  return(objMOList)
}


#' Performing peak annotation and motif enrichment analysis
#'
#' @aliases annotateATACPeaksMotif
#'
#' @description This function performs peak annotation and motif enrichment
#'              analysis on the ATACseq peaks. The peak annotation is performed
#'              using the ChIPseeker package, and the motif enrichment analysis
#'              is performed using the monaLisa package.
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
#' @param tssRegion A numeric vector of length 2, specifying the region around
#'                  the TSS to annotate with (default: +/- 3000). The first and
#'                  second elements of the vector specify the base pairs
#'                  upstream and downstream of the TSS, respectively.
#' @param TxDb The TxDb object to use for annotation, using one of the TxDb
#'             packages. See examples in how to generate the TxDb object.
#' @param annoDb The annotation database to use for annotation, using one of
#'               the org.*.eg.db packages. Must be a valid string for the name
#'               of the package. In addition, before calling the function, the
#'               user must load the corresponding annotation package using
#'              'library(org.*.eg.db)'. For example, if the user wants to use
#'              the org.Hs.eg.db package for human samples, the user must first
#'              load the package using 'library(org.Hs.eg.db)' and then supply
#'              the string "org.Hs.eg.db" to the annoDb argument. For details,
#'              see the vignette of the ChIPseeker package for
#'              custom annotation databases.
#' @param bsgenome The BSgenome object that contains the genome sequence.
#'                 Should be a valid BSgenome object generated using the
#'                 one of the BSgenome packages. The genome sequence should
#'                 correspond to the genome version used for the TxDb object.
#'                 Read the Bioconductor vignette for the BSgenome package for
#'                 more details.
#' @param pwmL A PWMatrixList object containing the position-weight matrices
#'             (PWMs) to use for motif enrichment analysis. The PWMs can be
#'             generated using the getMatrixSet() function from the TFBSTools
#'             package. See the vignette of the TFBSTools package for more
#'             details.
#' @param fixedWidth The fixed width to adjust the peak size for motif if
#'                   necessary (default: 500)
#'
#' @return An object of class MOList, with the ATACseq peaks annotated with
#'         genomic features and motif enrichment analysis results appended.
#'
#' @importFrom ChIPseeker annotatePeak
#' @importFrom GenomicRanges makeGRangesFromDataFrame resize
#' @importFrom monaLisa calcBinnedMotifEnrR
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate position-weight matrices (PWMs) from the JASPAR database
#' # This is an example of how to generate the PWMs from the JASPAR database
#' # Alternative databases can be used, as long as a PWMatrixList object is
#' # generated
#' # This example requires users to have the JASPAR2022 package installed,
#' # which is only available on Bioconductor 3.18 or later
#' # To install the package, run the following command in R:
#' # BiocManager::install("JASPAR2022")
#'
#' # Generate the PWMs
#' pwmL <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022,
#'   opts = list(
#'     matrixtype = "PWM",
#'     tax_group = "vertebrates"
#'   )
#' )
#' # Alternatively, the package has provided a pre-generated PWM list
#' data("jasparVertebratePWM")
#'
#' # Make use of the package's provided example data
#' data("expMOList")
#'
#' # Annotate the ATACseq peaks with genomic features and perform motif
#' # enrichment analysis
#' expMOList <- annotateATACPeaksMotif(expMOList,
#'   tssRegion = c(-3000, 3000),
#'   TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
#'   annoDb = "org.Hs.eg.db",
#'   bsgenome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
#'   pwmL = jasparVertebratePWM,
#'   fixedWidth = 500
#' )
#' }
#'
annotateATACPeaksMotif <- function(objMOList,
                                   tssRegion = c(-3000, 3000),
                                   TxDb,
                                   annoDb = "org.Hs.eg.db",
                                   bsgenome,
                                   pwmL,
                                   fixedWidth = 500) {
  # Annotate the ATACseq peaks with genomic features
  objMOList <- annotatePeaks(objMOList,
    tssRegion = tssRegion,
    TxDb = TxDb,
    annoDb = annoDb
  )

  # Perform motif enrichment analysis on the ATACseq peaks
  objMOList <- enrichMotifs(objMOList,
    bsgenome = bsgenome,
    pwmL = pwmL,
    fixedWidth = fixedWidth
  )
  return(objMOList)
}


#' Select enriched motifs
#'
#' @keywords internal
#'
#' @description This function selects enriched motifs based on user-specified
#'              cutoff and the motif enrichment result. If pValue is specified,
#'              pValueAdj will be ignored.
#'
#' @param enrichedMotifs A SummarizedExperiment object containing the motif
#'                       enrichment results
#' @param pValueAdj The cutoff for adjusted p-value, default is 0.05
#' @param pValue The cutoff for p-value, default is NULL. if pValue is
#'               specified, pValueAdj will be ignored
#' @param log2FEnrich The cutoff for log2 fold enrichment, default is NULL
#'
#' @return A logical vector with same number of rows as the number of motifs
#'         in the input object
#'
#' @importFrom SummarizedExperiment assay
#'
#' @references
#' \insertRef{SummarizedExperiment}{IntegraTRN}
#'
selectedMotifs <- function(enrichedMotifs,
                           pValueAdj = 0.05,
                           pValue = NULL,
                           log2FEnrich = NULL) {
  # Select by P-value or adjusted P-value
  if (is.null(pValue)) {
    # Use adjusted p-value
    p <- -log(pValueAdj, 10)
    selected <- apply(
      SummarizedExperiment::assay(
        enrichedMotifs,
        "negLog10Padj"
      ),
      1,
      function(x) max(abs(x), 0, na.rm = TRUE)
    ) > p
  } else {
    # Use p-value
    p <- -log(pValue, 10)
    selected <- apply(
      SummarizedExperiment::assay(
        enrichedMotifs,
        "negLog10P"
      ),
      1,
      function(x) max(abs(x), 0, na.rm = TRUE)
    ) > p
  }
  # Select by log2 fold enrichment
  if (!is.null(log2FEnrich)) {
    selected <- selected & apply(
      SummarizedExperiment::assay(
        enrichedMotifs,
        "log2enr"
      ),
      1,
      function(x) max(abs(x), 0, na.rm = TRUE)
    ) >
      log2FEnrich
  } else {
    # Do nothing
  }
  return(selected)
}


# [END]
