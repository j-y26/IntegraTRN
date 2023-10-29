# Purpose: Creating the Multi-Omics List S4 class
# Author: Jielin Yang
# Date: 2023-10-29
# Version: 1.0
# Bugs and Issues: None

# Define the MOList S4 class
# Slots: RNAseq, RNAseqSamples, smallRNAseq, smallRNAseqSamples, proteomics,
#        proteomicsSamples, ATACpeaks
# Inheritance: list
setClass("MOList",
  # Inheritance
  contains = "list",
  # Slots
  slots = list(
    RNAseq = "numeric",
    RNAseqSamples = "list",
    smallRNAseq = "numeric",
    smallRNAseqSamples = "list",
    proteomics = "numeric",
    proteomicsSamples = "list",
    ATACpeaks = "list"
  ),
  # Prototypes
  prototype = prototype(
    RNAseq = NA_real_,
    RNAseqSamples = list(
      samples = NA,
      groupBy = NA
    ),
    smallRNAseq = NA_real_,
    smallRNAseqSamples = list(
      samples = NA,
      groupBy = NA
    ),
    proteomics = NA_real_,
    proteomicsSamples = list(
      samples = NA,
      groupBy = NA
    ),
    ATACpeaks = list(
      peaksCond1 = NA,
      peaksCond2 = NA
    )
  )
)








# [END]
