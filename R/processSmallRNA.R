# Purpose: Small RNA processing functions, annotations, and clustering
# Author: Jielin Yang
# Date: 2023-11-01
# Version: 1.0
# Bugs and Issues: None


# Define global variables
HUMAN <- "human"
USERANNO <- "userAnno"
SMALLRNA_CATEGORIES <- c("miRNA", "piRNA", "tRNA", "circRNA", "snRNA", "snoRNA")

#' Check annotation coverage
#'
#' @keywords internal
#'
#' @description This function checks if the annotation covers all transcript
#'              names provided
#'
#' @param transcripts A character vector of transcript names
#' @param objMOList A object of class MOList
#' @param anno A character string indicating the annotation to be used. Either
#'            "human" or "userAnno"
#'
#' @return A boolean indicating whether the annotation covers all transcript
#'
#' @importFrom dplyr %>%
#'
#' @references
#' \insertRef{dplyr}{IntegraTRN}
#'
sncAnnoCoverage <- function(transcripts, objMOList, anno) {
  # Check if the annotation covers all transcripts
  if (anno == HUMAN) {
    passed <- transcripts %in% unlist(SNCANNOLIST_HSAPIENS)
  } else if (anno == USERANNO) {
    passed <- transcripts %in% unlist(objMOList$annoSncRNA)
  } else {
    stop("Invalid annotation type. Please use either \"human\" or
    \"userAnno\".")
  }

  return(all(passed))
}


#' Check if the annotation covers all small RNA transcripts in the MOList
#' object
#'
#' @keywords internal
#'
#' @description This function checks if the annotation covers all small RNA
#'              transcripts in the MOList object. If not, it will throw an
#'              warning message, but the function will still continue to run.
#'
#' @keywords internal
#'
#' @param objMOList A object of class MOList
#' @param anno A character string indicating the annotation to be used. Either
#'             "human" or "userAnno".
#'
#' @return NULL
#'
#' @importFrom dplyr %>% filter pull
#'
#' @references
#' \insertRef{dplyr}{IntegraTRN}
#'
checkSmallAnnoCoverage <- function(objMOList, anno) {
  # Check if small RNA data is available
  if (is.null(getRawData(objMOList, SMALLRNA))) {
    return(invisible(NULL))
  }
  # Otherwise, perform the check
  # Retrieve small RNA transcripts
  smallRNA <- rownames(getRawData(objMOList, SMALLRNA))

  # Check if the annotation covers all small RNA transcripts
  passed <- sncAnnoCoverage(smallRNA, objMOList, anno)

  # Throw warning message if not all small RNA transcripts are covered
  if (!passed) {
    warning("Not all small RNA transcripts are covered by the annotation,
             consider using a more comprehensive annotation.")
  } else {
    return(invisible(NULL))
  }

  # Check if differentially expressed smallRANs are covered, code executed
  # only if not all small RNA transcripts are covered
  if (is.null(objMOList$DEsmallRNAseq)) {
    warning("Differential expression analysis has not been performed. Consider
             re-running the analysis after performing differential expression
             analysis to check for annotation coverage on differentially
            expressed small RNAs.")
    return(invisible(NULL))
  }

  # Retrieve small RNA transcripts with adjusted p-value < 0.05
  # Code only executed if not all small RNA transcripts are covered
  deSmallRNA <- objMOList$DEsmallRNAseq %>%
    exportDE() %>%
    dplyr::filter(padj < 0.05) %>%
    rownames()
  dePassed <- sncAnnoCoverage(deSmallRNA, objMOList, anno)

  # Throw warning message if not all differentially expressed small RNA
  # transcripts are covered
  if (!dePassed) {
    warning("Not all small RNA transcripts with adjusted P-value < 0.05 are
             covered by the annotation, consider using a more comprehensive
             annotation.")
  } else {
    return(invisible(NULL))
  }
}


#' Extract small RNA categories
#'
#' @keywords internal
#'
#' @description This function extracts the small RNAs belonging to each of the
#'              categories, based on input.
#'
#' @param annoDF A data frame containing the annotation information for small
#'              RNA transcripts. The first column should be the transcript name
#'              of the small RNAs, and the second column should be the category
#'              of the small RNAs.
#' @param category A character string indicating the category of small RNAs to
#'                be extracted.
#'
#' @return A character vector containing the small RNA transcripts belonging to
#'        the specified category.
#'
#' @importFrom dplyr %>% filter pull
#'
#' @references
#' \insertRef{dplyr}{IntegraTRN}
#'
extractTranscriptFromAnno <- function(annoDF, category) {
  # Extract small RNA transcripts belonging to the specified category
  transcripts <- annoDF[grepl(category, annoDF[, 2], ignore.case = TRUE), 1]

  # Check if it is empty
  if (length(transcripts) == 0) {
    warning("No small RNA transcripts of ", category, " are found.")
    return(NULL)
  } else {
    return(transcripts)
  }
}


#' Annotate small RNA transcripts to respective categories
#'
#' @aliases annotateSmallRNA
#'
#' @description This function annotates small RNA transcripts to their
#'              respective categories. The categories are defined as follows:
#'              1. miRNA
#'              2. piRNA
#'              3. tRNA
#'              4. circRNA
#'              5. snRNA
#'              6. snoRNA
#'
#' @param objMOList A object of class MOList
#' @param anno A data frame containing
#'             annotation information for small RNA transcripts. The first
#'             column should be the transcript name of the small RNAs, and
#'             the second column should be the category of the small RNAs.
#'             The category should be one of the following: miRNA, piRNA,
#'             tRNA, circRNA, snRNA, snoRNA. Or, the input can be a simple
#'             string "human" as a shortcut to use the default annotation
#'             provided internally with the package.
#' \itemize{
#'  \item{transcript}{The transcript name of the small RNA}
#'  \item{category}{The category of the small RNA}
#' }
#'
#' @details The package provides an internal annotation for human small RNAs
#'          based on the hg38 genome. The annotation was generated by COMPSRA,
#'          using annotations from the following databases:
#'             1. miRbase
#'             2. piRNABank
#'             3. piRBase
#'             4. GtRNAdb
#'             5. circBase
#'             6. GENCODE
#'             7. piRNACluster
#'             Please see the reference for more details on how the annotation
#'             was generated.
#'
#'
#' @return A MOList object with the annotation information added to the
#'        annoSncRNA element.
#'
#' @references
#' \insertRef{li2020compsra}{IntegraTRN}
#'
#' \insertRef{kozomara2010mirbase}{IntegraTRN}
#'
#' \insertRef{chan2016gtrnadb}{IntegraTRN}
#'
#' \insertRef{harrow2012gencode}{IntegraTRN}
#'
#' \insertRef{glavzar2014circbase}{IntegraTRN}
#'
#' \insertRef{sai2008pirnabank}{IntegraTRN}
#'
#' \insertRef{zhang2014pirbase}{IntegraTRN}
#'
#' \insertRef{rosenkranz2016pirna}{IntegraTRN}
#'
#' \insertRef{dplyr}{IntegraTRN}
#'
#' @export
#'
#' @importFrom dplyr %>% filter pull
#'
#' @examples
#' # Example 1: Adding user-defined annotation
#' # Generate some sample RNAseq and smallRNAseq data
#' rnaMatrix <- matrix(sample(0:100, 1000, replace = TRUE), nrow = 100)
#' rownames(rnaMatrix) <- paste0("gene_", seq_len(100))
#' rnaGroupBy <- rep(c("group1", "group2"), each = 5)
#' smallRNAMatrix <- matrix(sample(0:100, 1000, replace = TRUE), nrow = 100)
#' rownames(smallRNAMatrix) <- paste0("transcript_", seq_len(100))
#' smallRNAGroupBy <- rep(c("group1", "group2"), each = 5)
#'
#' # Create example MOList object
#' objMOList <- MOList(
#'   RNAseq = rnaMatrix, RNAGroupBy = rnaGroupBy,
#'   smallRNAseq = smallRNAMatrix, smallRNAGroupBy = smallRNAGroupBy
#' )
#'
#' # Create example annotation
#' anno <- data.frame(
#'   transcript = paste0("transcript_", seq_len(120)),
#'   category = rep(c("miRNA", "piRNA", "tRNA", "circRNA", "snRNA", "snoRNA"),
#'     each = 20
#'   )
#' )
#'
#' # Annotate small RNA transcripts
#' annotateSmallRNA(objMOList, anno)
#'
#' # Example 2: Using the preprocessed annotation provided by the package
#' # This time we use the package provided example data
#' data(expMOList)
#'
#' # Annotate small RNA transcripts
#' annotateSmallRNA(expMOList, anno = "human")
#'
#' # or simply use the default parameter
#' annotateSmallRNA(expMOList)
#'
annotateSmallRNA <- function(objMOList, anno = "human") {
  # Construct annotation based on the input
  if (all(anno == HUMAN)) {
    # Use the internal annotation
    objMOList$annoSncRNA <- HUMAN
  } else {
    userAnno <- anno

    # Check of valid input as a data frame or a string
    if (inherits(userAnno, "data.frame")) {
      # Do nothing
    } else {
      stop("Invalid annotation input. Please provide a valid path to a
      tab-delimited file or a loaded data frame. See ?annotateSmallRNA for
      more details.")
    }

    # Generate annotation
    miRNAanno <- extractTranscriptFromAnno(userAnno, "miRNA")
    piRNAanno <- extractTranscriptFromAnno(userAnno, "piRNA")
    tRNAanno <- extractTranscriptFromAnno(userAnno, "tRNA")
    circRNAanno <- extractTranscriptFromAnno(userAnno, "circRNA")
    snRNAanno <- extractTranscriptFromAnno(userAnno, "snRNA")
    snoRNAanno <- extractTranscriptFromAnno(userAnno, "snoRNA")
    annoSncRNA <- list(
      miRNA = miRNAanno, piRNA = piRNAanno,
      tRNA = tRNAanno, circRNA = circRNAanno,
      snRNA = snRNAanno, snoRNA = snoRNAanno
    )
    objMOList$annoSncRNA <- annoSncRNA
    anno <- USERANNO
  }

  # Check if the annotation covers all small RNA transcripts

  checkSmallAnnoCoverage(objMOList, anno)

  # Return the object
  return(objMOList)
}


# [END]
