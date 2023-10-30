#' Gene expression of fetal heart tissues
#'
#' An RNAseq experiment conducted on fetal ventricular heart tissues from 53
#' individuals
#'
#' @format A matrix with 13043 rows (genes) and 53 columns (samples)
#' \describe{
#'  \item{sample_1, sample_2, sample_3, sample_4, sample_5, sample_6, sample_7,
#'  sample_8, sample_9, sample_10, sample_11, sample_12, sample_13, sample_14,
#'  sample_15, sample_16, sample_17, sample_18, sample_19, sample_20, sample_21,
#'  sample_22, sample_23, sample_24, sample_25, sample_26, sample_27, sample_28,
#'  sample_29, sample_30, sample_31, sample_32, sample_33, sample_34, sample_35,
#'  sample_36, sample_37, sample_38, sample_39, sample_40, sample_41, sample_42,
#'  sample_43, sample_44, sample_45, sample_46, sample_47, sample_48, sample_49,
#'  sample_50, sample_51, sample_52, sample_53}{Gene expression values for each
#'  sample}
#' }
#'
#' @source Paul Delgado Olguin Lab, The Hospital for Sick Children, Toronto,
#' ON, Canada
#'
#' @examples
#' \dontrun{
#' RNAseq_heart
#' }
#'
#' @references
#' Ahmed, A., Liang, M., Chi, L., Zhou, Y. Q., Sled, J. G., Wilson, M. D., &
#' Delgado-Olguín, P. (2021). Maternal obesity persistently alters cardiac
#' progenitor gene expression and programs adult-onset heart disease
#' susceptibility. Molecular Metabolism, 43.
#' https://doi.org/10.1016/j.molmet.2020.101116
"RNAseq_heart"

#' Sample information for RNAseq of heart tissues
#'
#' A selected set of sample information for the RNA-seq experiment conducted on
#' fetal ventricular heart tissues from 53 individuals
#'
#' @format A data frame with 53 rows (samples) and 4 columns (fields)
#' \describe{
#'  \item{Sample}{Sample name corresponding to the RNAseq gene expression
#'  dataset}
#'  \item{Age}{Age of the fetus in gestational weeks}
#'  \item{Sex}{Sex of the fetus}
#'  \item{Batch}{Batch number of the RNAseq experiment, represented from
#'  A to F}
#' }
#'
#' @source Paul Delgado Olguin Lab, The Hospital for Sick Children, Toronto,
#' ON, Canada
#'
#' @examples
#' \dontrun{
#' RNAseq_heart_samples
#' }
#'
#' @references
#' Ahmed, A., Liang, M., Chi, L., Zhou, Y. Q., Sled, J. G., Wilson, M. D., &
#' Delgado-Olguín, P. (2021). Maternal obesity persistently alters cardiac
#' progenitor gene expression and programs adult-onset heart disease
#' susceptibility. Molecular Metabolism, 43.
#' https://doi.org/10.1016/j.molmet.2020.101116
"RNAseq_heart_samples"

#' Expression of small RNA in fetal heart tissues
#'
#' A small RNAseq experiment conducted on fetal ventricular heart tissues from
#' 37 individuals
#'
#' @format A data frame with 14801 rows (small RNA transcripts) and 37 columns
#' (samples)
#' \describe{
#'  \item{sample_1, sample_2, sample_3, sample_4, sample_5, sample_6, sample_7,
#'  sample_8, sample_9, sample_10, sample_11, sample_12, sample_13, sample_14,
#'  sample_15, sample_16, sample_17, sample_18, sample_19, sample_20, sample_21,
#'  sample_22, sample_23, sample_24, sample_25, sample_26, sample_27, sample_28,
#'  sample_29, sample_30, sample_31, sample_32, sample_33, sample_34, sample_35,
#'  sample_36, sample_37}{Transcript expression values for each sample}
#' }
#'
#' @source Paul Delgado Olguin Lab, The Hospital for Sick Children, Toronto,
#' ON, Canada
#'
#' @examples
#' \dontrun{
#' smallRNA_heart
#' }
"smallRNAseq_heart"

#' Sample information for small RNAseq of heart tissues
#'
#' A selected set of sample information for the small RNA-seq experiment
#' conducted on fetal ventricular heart tissues from 37 individuals
#'
#' @format A data frame with 37 rows (samples) and 3 columns (fields)
#' \describe{
#'  \item{Sample}{Sample name corresponding to the small RNAseq gene expression
#'  dataset}
#'  \item{Age}{Age of the fetus in gestational weeks}
#'  \item{Sex}{Sex of the fetus}
#' }
#'
#' @source Paul Delgado Olguin Lab, The Hospital for Sick Children, Toronto,
#' ON, Canada
#'
#' @examples
#' \dontrun{
#' smallRNA_heart_samples
#' }
"smallRNAseq_heart_samples"

#' Protein abundance in fetal heart tissues
#'
#' A proteomics experiment conducted on fetal ventricular heart tissues from
#' 6 individuals, using LC-MS/MS. Only concensus proteins found in all samples
#' in the same conditions are included.
#'
#' @format A data frame with 489 rows (proteins) and 6 columns (samples)
#' \describe{
#'  \item{sample_1, sample_2, sample_3, sample_4, sample_5, sample_6}{Protein
#'  abundance values for each sample}
#' }
#'
#' @source Paul Delgado Olguin Lab, The Hospital for Sick Children, Toronto,
#' ON, Canada
#'
#' @examples
#' \dontrun{
#' protein_heart
#' }
"protein_heart"

#' Sample information for proteomics of heart tissues
#'
#' A selected set of sample information for the proteomics experiment conducted
#' on fetal ventricular heart tissues from 6 individuals
#'
#' @format A data frame with 6 rows (samples) and 3 columns (fields)
#' \describe{
#'  \item{Sample}{Sample name corresponding to the proteomics gene expression
#'  dataset}
#'  \item{Age}{Age of the fetus in gestational weeks}
#' }
#'
#' @source Paul Delgado Olguin Lab, The Hospital for Sick Children, Toronto,
#' ON, Canada
#'
#' @examples
#' \dontrun{
#' protein_heart_samples
#' }
"protein_heart_samples"
