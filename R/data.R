#' Gene expression of fetal heart tissues
#'
#' An RNAseq experiment conducted on fetal ventricular heart tissues from 53
#' individuals
#'
#' @format A matrix with 5000 rows (genes) and 53 columns (samples)
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
#' \insertRef{ahmed2021maternal}{IntegraTRN}
#'
"RNAseq_heart"

#' Sample information for RNAseq of heart tissues
#'
#' A selected set of sample information for the RNA-seq experiment conducted on
#' fetal ventricular heart tissues from 53 individuals
#'
#' @format A data frame with 53 rows (samples) and 4 columns (fields)
#' \describe{
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
#' \insertRef{ahmed2021maternal}{IntegraTRN}
#'
"RNAseq_heart_samples"

#' Expression of small RNA in fetal heart tissues
#'
#' A small RNAseq experiment conducted on fetal ventricular heart tissues from
#' 37 individuals
#'
#' @format A data frame with 1522 rows (small RNA transcripts) and 37 columns
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
#'
#' @references
#' \insertRef{adar2023integrate}{IntegraTRN}
#'
"smallRNAseq_heart"

#' Sample information for small RNAseq of heart tissues
#'
#' A selected set of sample information for the small RNA-seq experiment
#' conducted on fetal ventricular heart tissues from 37 individuals
#'
#' @format A data frame with 37 rows (samples) and 3 columns (fields)
#' \describe{
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
#'
#' @references
#' \insertRef{adar2023integrate}{IntegraTRN}
#'
"smallRNAseq_heart_samples"

#' Protein abundance in fetal heart tissues
#'
#' A proteomics experiment conducted on fetal ventricular heart tissues from
#' 6 individuals, using LC-MS/MS. Only concensus proteins found in all samples
#' in the same conditions are included. Some simulated values are included to
#' increase protein coverage to demonstrate package functionality.
#'
#' @format A data frame with 1448 rows (proteins) and 6 columns (samples)
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
#'
#' @references
#' \insertRef{adar2023integrate}{IntegraTRN}
#'
"protein_heart"

#' Sample information for proteomics of heart tissues
#'
#' A selected set of sample information for the proteomics experiment conducted
#' on fetal ventricular heart tissues from 6 individuals
#'
#' @format A data frame with 6 rows (samples) and 3 columns (fields)
#' \describe{
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
#'
#' @references
#' \insertRef{adar2023integrate}{IntegraTRN}
#'
"protein_heart_samples"

#' Example miRNA - mRNA interactions
#'
#' A list of miRNA-target interactions on fetal heart differentially expressed
#' genes. The list contains 2 vectors of the same length. Each representing a
#' pair of miRNA-target interactions. Curated partially from randomly selected
#' miRNA-target interactions from miRNet 2.0. Partially simulated.
#'
#' @format A list of 2 vectors of the same length. Each representing a pair of
#'         miRNA-target interactions.
#' \describe{\
#' \item{regulator}{A vector of miRNAs}
#' \item{target}{A vector of upregulated fetal heart mRNAs}
#' }
#'
#'
#' @references
#' \insertRef{chang2020mirnet}{IntegraTRN}
#'
#' @examples
#' \dontrun{
#' miR2Genes
#' }
#'
"miR2Genes"

#' Example TF - mRNA interactions
#'
#' A list of TF-target interactions on fetal heart differentially expressed
#' genes. The list contains 2 vectors of the same length. Each representing a
#' pair of TF-target interactions. Curated partially from randomly selected
#' TF-target interactions from TRRUST v2 using miRNet 2.0. Partially simulated
#' using simulated ATACseq peaks.
#'
#' @format A list of 2 vectors of the same length. Each representing a pair of
#'        TF-target interactions.
#' \describe{\
#' \item{regulator}{A vector of TFs}
#' \item{target}{A vector of upregulated fetal heart mRNAs}
#' }
#'
#' @references
#' \insertRef{chang2020mirnet}{IntegraTRN}
#'
#' @examples
#' \dontrun{
#' tf2Genes
#' }
#'
"tf2Genes"

#' JASPAR curated position weight matrices of transcription factor binding sites
#' for vertebrates based on database JASPAR CORE 2022
#'
#' A PWMatrixList object containing the position weight matrices of
#' transcription factor binding sites for vertebrates based on database JASPAR
#' CORE 2022. This list is extracted from the JASPAR2022 package using the
#' function \code{getMatrixSet()} from the \code{TFBSTools} package.
#'
#' @format A list of 841 position weight matrices of transcription factor
#' binding sites for vertebrates, as a PWMatrixList object
#'
#' @source JASPAR2022 package
#'
#' @references
#' \insertRef{castro2022jaspar}{IntegraTRN}
#'
#' \insertRef{tan2016tfbstools}{IntegraTRN}
#'
#' @examples
#' \dontrun{
#' jasparVertebratePWM
#' }
#'
"jasparVertebratePWM"

#' Protein to gene mapping
#'
#' A data frame containing the mapping between protein names and gene names
#' used in the example proteomics and RNAseq datasets. This is only a subset
#' of the full mapping and does not cover all proteins in the example data for
#' illustration purposes.
#'
#' @format A data frame with 1421 rows (proteins) and 2 columns (fields)
#' \describe{
#' \item{protein}{Protein names}
#' \item{gene}{Gene names}
#' }
#'
#' @references
#' \insertRef{durinck2009mapping}{IntegraTRN}
#'
#' @examples
#' \dontrun{
#' proteinGeneIDConvert
#' }
#'
"proteinGeneIDConvert"

#' Example MOList S4 object
#'
#' An example MOList S4 object containing the RNAseq, small RNAseq, proteomics,
#' and ATACseq peaks
#'
#' @format An S4 object of class MOList
#' \describe{
#' \item{RNAseq}{A matrix of RNAseq gene expression values with 100 genes}
#' \item{smallRNAseq}{A matrix of small RNAseq gene expression values with 100
#'                    genes}
#' \item{proteomics}{A matrix of proteomics gene expression values with 100
#'                   genes}
#' \item{ATACseq}{A list of two data frames of ATACseq peaks for 2 conditions,
#'               each with 100 peaks}
#' \item{DERNAseq}{A DETag object of the differentially expressed genes from
#'                 RNAseq}
#' \item{DEsmallRNAseq}{A DETag object of the differentially expressed genes
#'                     from small RNAseq}
#' \item{DEproteomics}{A DETag object of the differentially expressed genes
#'                    from proteomics}
#' \item{DEATAC}{A DETag object of peaks}
#' \item{anncSncRNA}{The string "human"}
#' \item{matchingRNAsmallRNA}{A list of matched indices between RNAseq and
#'                           small RNAseq samples}
#' \item{extInteractions}{A list of 2 lists for miRNA-target and TF-target
#'                       interactions. Each contains 50 pairs of interactions}
#' }
#' @source Paul Delgado Olguin Lab, The Hospital for Sick Children, Toronto,
#' ON, Canada
#'
#' @references
#' \insertRef{ahmed2021maternal}{IntegraTRN}
#'
#' \insertRef{adar2023integrate}{IntegraTRN}
#'
#' @examples
#' \dontrun{
#' expMOList
#' }
#'
"expMOList"

#' Small RNA type annotation for human
#'
#' A list of 6 character vectors of small RNA. Each element of the list
#' represents a small RNA type. The 6 types are: miRNA, piRNA, snoRNA, snRNA,
#' tRNA, and circRNA. The character vectors contain the names of small RNAs
#' belonging to the corresponding type, using the nomenclature from the
#' listed databases.
#'
#' @format A list of 6 character vectors of small RNA
#' \describe{
#' \item{miRNA}{A character vector of miRNA names}
#' \item{piRNA}{A character vector of piRNA names}
#' \item{snoRNA}{A character vector of snoRNA names}
#' \item{snRNA}{A character vector of snRNA names}
#' \item{tRNA}{A character vector of tRNA names}
#' \item{circRNA}{A character vector of circRNA names}
#' }
#'
#' @source miRbase, piRNAbank, piRBase, GtRNAdb, circBase, GENCODE, and
#'         piRNACluster
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
#' @examples
#' \dontrun{
#' sncRNAAnnotation
#' }
"sncRNAAnnotation"

# [END]
