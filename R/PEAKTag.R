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
#' \item \code{Condition}: Condition of the peak, "-" for peaks from condition 1
#'                         and "+" for peaks from condition 2
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
#' @references
#' \insertRef{yu2015chipseeker}{IntegraTRN}
#'
setClass("PEAKTag",
  representation(
    annotatedPeaks = "csAnno",
    TxDB = "TxDb",
    annoDB = "character"
  ),
  contains = c("DETag", "list")
)


#' @title Constructor for the PEAKTag class
#'
#' @aliases PEAKTag
#'
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
#' @importClassesFrom ChIPseeker csAnno
#'
#' @export
#'
#' @references
#' \insertRef{yu2015chipseeker}{IntegraTRN}
#'
#' \insertRef{TxDb.Hsapiens.UCSC.hg38.knownGene}{IntegraTRN}
#'
#' \insertRef{org.Hs.eg.db}{IntegraTRN}
#'
#' @examples
#' # Make use of the package's provided example data
#' data(expMOList)
#'
#' # Retrieve the DETag object for the ATACseq peaks
#' atacDETag <- expMOList$DEATAC
#'
#' # Example 1: Construct a PEAKTag object without annotation
#' atacPEAKTag <- PEAKTag(atacDETag)
#'
#' # Example 2: Construct a PEAKTag object with annotation
#' # Suppose TxDB and annoDB are already loaded
#' \dontrun{
#' TxDB <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#' atacPEAKTag <- PEAKTag(atacDETag, TxDB = TxDB, annoDB = "org.Hs.eg.db")
#' }
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
  objPEAKTag <- as(objDETag, "PEAKTag")
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

#' print PEAKTag
#'
#' @aliases show,PEAKTag-method
#'
#' @description The print method for the PEAKTag class.
#'
#' @param object A PEAKTag object
#'
#' @return NULL
#'
#' @export
#'
#' @examples
#' # Use the package's provided example data
#' data(expMOList)
#'
#' # Print the object
#' print(expMOList$DEATAC)
#'
#' # Or simply type the object name
#' expMOList$DEATAC
#'
setMethod(
  "show", "PEAKTag",
  function(object) {
    cat(
      "A PEAKTag object containing", nrow(object@DEResult),
      "differential peaks.\n"
    )
    if (length(object@annoDB) > 0) {
      cat("Annotation databases: ", object@annoDB, "\n")
      cat("Transcript databases: ", object@TxDB$packageName, "\n\n")
    } else {
      cat("No available peak annotation.\n")
    }
    cat(
      "Number of peaks in condition 1: ",
      sum(object@DEResult$Condition == "-"), "\n"
    )
    cat(
      "Number of peaks in condition 2: ",
      sum(object@DEResult$Condition == "+"), "\n\n"
    )
    cat("A snapshot of the peaks:\n")
    as.data.frame(object) %>%
      head(5) %>%
      print()
    return(invisible(NULL))
  }
)


#' Convert to a GRanges object
#'
#' @aliases asGRanges,PEAKTag-method
#'
#' @export
#'
#' @importFrom ChIPseeker as.GRanges
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#'
#' @param x A PEAKTag object
#'
#' @return A GRanges object
#'
#' @examples
#' # Use the package's provided example data
#' data(expMOList)
#'
#' # Convert the object to a GRanges object
#' # Requires that annotation has been performed
#' \dontrun{
#' asGRanges(expMOList$DEATAC)
#' }
#'
setGeneric(
  "asGRanges",
  function(x) {
    standardGeneric("asGRanges")
  }
)
setMethod(
  "asGRanges",
  signature(x = "PEAKTag"),
  function(x) {
    # Determine if an annotation exists
    if (length(x@annotatedPeaks@hasGenomicAnnotation) == 0) {
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
#' @aliases as.data.frame,PEAKTag-method
#'
#' @export
#'
#' @importFrom BiocGenerics as.data.frame
#'
#' @param x A PEAKTag object
#'
#' @return A data frame
#'
#' @examples
#' # Use the package's provided example data
#' data(expMOList)
#'
#' # Convert the object to a data frame
#' # Requires that annotation has been performed
#' \dontrun{
#' as.data.frame(expMOList$DEATAC)
#' }
#'
setMethod(
  "as.data.frame",
  signature(x = "PEAKTag"),
  function(x) {
    # Determine if an annotation exists
    if (length(x@annotatedPeaks@hasGenomicAnnotation) == 0) {
      # No annotation
      peakDF <- x@DEResult
    } else {
      # Annotation exists
      peakDF <- as.data.frame(x@annotatedPeaks)
    }
    return(peakDF)
  }
)


# [END]
