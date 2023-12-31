% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diffOmicsAnalysis.R
\name{filterGeneCounts}
\alias{filterGeneCounts}
\title{Filter the gene counts for the count-based omics data}
\usage{
filterGeneCounts(objMOList, omic)
}
\arguments{
\item{objMOList}{A MOList object containing the omics data}

\item{omic}{A character string specifying the omics data to be filtered
must be one of "RNAseq", "smallRNAseq", and "proteomics"}
}
\value{
An MOList object containing the filtered omics data
\itemize{
\item code(RNAseq): The RNAseq data in the MOList object is filtered and
                    updated, if omic = "RNAseq"
\item code(smallRNAseq): The small RNAseq data in the MOList object is
                          filtered and updated, if omic = "smallRNAseq"
\item code(proteomics): The protein data in the MOList object is filtered and
                        updated, if omic = "proteomics"
\item code{ATACpeaks}: A list containing the ATAC peaks for condition 1 and
                       condition 2, kept unchanged
}
}
\description{
Filtering is based on the design of the experiment. If the
             samples are only grouped into 2 conditions, then the genes with
             counts in more than 1 counts per million (CPM) in at least
             min(#samples in condition 1, #samples in condition 2) samples
             are kept. If the samples are grouped by a continuous variable,
             then the genes with counts in more than 1 CPM in at least
             30% of the samples are kept.
}
\references{
\insertRef{robinson2010edger}{IntegraTRN}
}
\keyword{internal}
