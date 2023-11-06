# Purpose: A simple class to store the peak information
# Author: Jielin Yang
# Date: 2023-11-06
# Version: 1.0
# Bugs and Issues: None


#' @title PEAKTag S4 class for storing peak information
#'
#' @description The PEAKTag S4 class inherits the DETag class and is a more
#'              specialized class for storing peak information. It contains
#'              all slots of the DETag class, with additional slots added as a
#'              more convenient way for users to access the peak information.
#'
#' @exportClass PEAKTag
#' @importClassesFrom ChIPseeker csAnno
#'
#' @slot DEResult A data frame containing the annotated peak information
#' \itemize{
#'  \item \code{chr}: Chromosome
#' \item \code{start}: Start coordinate
#' \item \code{end}: End coordinate
#' \item \code{width}: Width of the peak
#' \item \code{strand}: Strand, always \code{"*"} for ATACseq
#' \item \code{Condition}: Condition of the peak, "-" for peaks from condition 1
#'                         and "+" for peaks from condition 2
#' \item \code{annotation}: Annotation of the peak, the class of the peak
#' \item \code{ENSEMBL}: ENSEMBL gene ID, based on the peak annotation
#' \item \code{SYMBOL}: Gene symbol, based on the peak annotation
#' \item \code{GENENAME}: Gene name, based on the peak annotation
#' }
#'
#' @slot method A character string indicating the method used to generate the
#'              peak set, defined similarly to the DETag class. In this case,
#'              must be "GRangesATAC"
#' @slot normalizedCounts A matrix containing the normalized counts for the
#'                        genes. Must be NULL for this class. However, this
#'                        slot is kept for compatibility with the DETag class
#'                        and allows future expansion of the package to support
#'                        differential peak calling using sequence coverage
#'                        information.
#' @slot annotatedPeaks A csAnno object containing the annotated peaks, as
#'                      defined in the ChIPseeker package
#' @slot TxDB A TxDb object containing the transcript database used for the
#'            annotation
#' @slot annoDB A string indicating the annotation database used for the
#'              annotation
#' @slot .Data A list containing additional peak information
#'
methods::setClass("PEAKTag",
  representation(
    annotatedPeaks = "csAnno",
    TxDB = "TxDb",
    annoDB = "character"
  ),
  contains = c("DETag", "list")
)


#' @title Constructor for the PEAKTag class
#' @description This function constructs a PEAKTag object from a DETag object
#'              and the annotated peaks, as defined in the ChIPseeker package.
#'
#' @param objDETag A DETag object
#' @param annotatedPeaks A csAnno object containing the annotated peaks, as
#'                       defined in the ChIPseeker package
#' @param TxDB A TxDb object containing the transcript database used for the
#'             annotation
#' @param annoDB A string indicating the annotation database used for the
#'               annotation
#'
#' @return A PEAKTag object
#'
#' @importClassesFrom DETag DETag
#' @importClassesFrom ChIPseeker csAnno
#' @importFrom methods as
#'
#' @export
#'
PEAKTag <- function(objDETag,
                    annotatedPeaks = NULL,
                    TxDB = NULL,
                    annoDB = NULL) {
  # Validating inputs
  if (!is(objDETag, "DETag")) {
    stop("The input must be a DETag object")
  }
  if (!is.null(annotatedPeaks) && !is(annotatedPeaks, "csAnno")) {
    stop("The annotated peaks must be a csAnno object")
  }
  if (!is.null(TxDB) && !is(TxDB, "TxDb")) {
    stop("The TxDb object must be a TxDb object")
  }
  if (!is.null(annoDB) && !is.character(annoDB)) {
    stop("The annotation database must be a character string")
  }

  # The DETag object must be generated using the GRangesATAC method
  if (objDETag@method != ATAC_GRANGE) {
    stop("Not an ATACseq peak set. Cannot coerce to PEAKTag object.")
  } else {
    # Continue
  }

  # Coerce the DETag object to a PEAKTag object
  objPEAKTag <- methods::as(objDETag, "PEAKTag")
  if (!is.null(annotatedPeaks)) {
    objPEAKTag@annotatedPeaks <- annotatedPeaks
  }
  if (!is.null(TxDB)) {
    objPEAKTag@TxDB <- TxDB
  }
  if (!is.null(annoDB)) {
    objPEAKTag@annoDB <- annoDB
  }

  return(objPEAKTag)
}


# Define some S4 methods for the PEAKTag class

#' @rdname PEAKTag-class
#'
#' @export
#'
#' @importFrom methods setMethod
#'
#' @method print PEAKTag
#'
#' @description The print method for the PEAKTag class.
#'
#' @inheritParams DETag-class
#'
#' @return NULL
#'
#' @examples
#' # Assuming that the object "peakTag" is a PEAKTag object
#' \dontrun{
#' print(peakTag)
#' }
#'
#' # Or simply type the object name
#' \dontrun{
#' peakTag
#' }
#'
methods::setMethod(
  "show", "PEAKTag",
  function(object) {
    cat(
      "A PEAKTag object containing", nrow(object@DEResult),
      "differential peaks.\n"
    )
    cat("Annotation databases: ", object@annoDB, "\n")
    cat("Transcript databases: ", object@TxDB$packageName, "\n\n")
    cat(
      "Number of peaks in condition 1: ",
      sum(object@DEResult$Condition == "-"), "\n"
    )
    cat(
      "Number of peaks in condition 2: ",
      sum(object@DEResult$Condition == "+"), "\n\n"
    )
    cat("A snapshot of the peaks:\n")
    anno <- object@annotatedPeaks
    if (is.null(anno)) {
      print(head(object@DEResult, 5))
    } else {
      peakDF <- csAnnoToDF(object@annotatedPeaks)
      print(head(peakDF, 5))
    }
    rm(anno, peakDF)
    return(invisible(NULL))
  }
)


#' Convert to a GRanges object
#'
#' @rdname PEAKTag-class
#'
#' @export
#'
#' @importFrom methods setMethod setGeneric
#' @importFrom ChIPseeker as.GRanges
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#'
#' @param x A PEAKTag object
#'
#' @return A GRanges object
#'
#' @examples
#' # Assuming that the object "peakTag" is a PEAKTag object
#' \dontrun{
#' asGRanges(peakTag)
#' }
#'
methods::setGeneric(
  "asGRanges",
  function(x) {
    standardGeneric("asGRanges")
  }
)
methods::setMethod(
  "asGRanges",
  signature(x = "PEAKTag"),
  function(x) {
    # Determine if an annotation exists
    if (is.null(x@annotatedPeaks)) {
      # No annotation
      peakGR <- GenomicRanges::makeGRangesFromDataFrame(
        x@DEResult,
        keep.extra.columns = TRUE
      )
    } else {
      # Annotation exists
      peakGR <- ChIPseeker::as.GRanges(x@annotatedPeaks)
    }
    return(peakGR)
  }
)


#' Convert to a data frame
#'
#' @rdname PEAKTag-class
#'
#' @export
#'
#' @importFrom methods setMethod setGeneric
#'





# [END]
