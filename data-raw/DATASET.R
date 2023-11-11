library(tidyverse)


# Internal data to be used

# Preparing human transcription factor data

humanTF <- read.csv("./data-raw/DatabaseExtract_v_1.01.csv",
  header = TRUE,
  row.names = 1
)
humanTF <- humanTF %>% dplyr::select(
  HGNC.symbol, Ensembl.ID, EntrezGene.ID,
  EntrezGene.Description, DBD, PDB.ID, Is.TF.
)

rownames(humanTF) <- humanTF$HGNC.symbol

TFHsapiens <- humanTF


# Document human small RNAs
write.csv(miRNAHsapiens, "./data-raw/miRNAHsapiens_COMPSRA.csv",
  row.names = FALSE
)
write.csv(snoRNAHsapiens, "./data-raw/snoRNAHsapiens_COMPSRA.csv",
  row.names = FALSE
)
write.csv(snRNAHsapiens, "./data-raw/snRNAHsapiens_COMPSRA.csv",
  row.names = FALSE
)
write.csv(tRNAHsapiens, "./data-raw/tRNAHsapiens_COMPSRA.csv",
  row.names = FALSE
)
write.csv(piRNAHsapiens, "./data-raw/piRNAHsapiens_COMPSRA.csv",
  row.names = FALSE
)
write.csv(circRNAHsapiens, "./data-raw/circRNAHsapiens_COMPSRA.csv",
  row.names = FALSE
)


# Saving internal data
usethis::use_data(SNCANNOLIST_HSAPIENS, internal = TRUE, overwrite = TRUE)


# Exported Data

# Preprocess RNAseq data
RNAseq <- read.table("./data-raw/GSE241758_allcounts.txt",
  header = TRUE, sep = "\t"
)

# Rename the columns
sampleNames <- paste0("sample_", seq(1, ncol(RNAseq) - 1))
colName <- c("gene", sampleNames)
colnames(RNAseq) <- colName

# Filtering low count genes
cpm <- edgeR::cpm(RNAseq[, -1])
rownames(cpm) <- RNAseq$gene
keep <- rowSums(cpm > 1) == ncol(cpm)
RNAseq <- RNAseq[keep, ]

# Save the RNAseq count data
write.csv(RNAseq, "./data-raw/RNAseq_heart.csv", row.names = FALSE)
RNAseq_heart <- RNAseq
usethis::use_data(RNAseq_heart, overwrite = TRUE)



# Sample matrix for RNAseq data
RNAseqSamples <- readxl::read_excel("./data-raw/Additional_file6.xlsx")
RNAseqSamples <- RNAseqSamples[, c("Age", "Gender", "Batch")]
colnames(RNAseqSamples)[2] <- "Sex"
RNAseqSamples <- cbind(sampleNames, RNAseqSamples)
colnames(RNAseqSamples)[1] <- "Sample"
write.csv(RNAseqSamples, "./data-raw/RNAseq_heart_samples.csv",
  row.names = FALSE
)
RNAseq_heart_samples <- RNAseqSamples
usethis::use_data(RNAseq_heart_samples, overwrite = TRUE)



# Preprocess small RNAseq data
smallRNAseq <- read.table(
  "./data-raw/GSE241759_all_unique_counts_smallRNAseq.txt",
  header = TRUE, sep = "\t"
)

# Filtering out genes that has more than half of the samples with 0 counts
keep <- rowSums(smallRNAseq > 0) > ncol(smallRNAseq) / 2
smallRNAseq <- smallRNAseq[keep, ]
sampleNames <- paste0("sample_", seq(1, ncol(smallRNAseq)))
colnames(smallRNAseq) <- sampleNames
smallRNAseq <- rownames_to_column(smallRNAseq, var = "transcript")
write.csv(smallRNAseq, "./data-raw/smallRNAseq_heart.csv", row.names = FALSE)
smallRNAseq_heart <- smallRNAseq
usethis::use_data(smallRNAseq_heart, overwrite = TRUE)



# Sample matrix for small RNAseq data
smallRNAseqSamples <- readxl::read_excel("./data-raw/Additional_file1.xlsx")
smallRNAseqSamples <- smallRNAseqSamples[
  ,
  c("Sample ID", "Gestational Age (wks)", "Sex")
]
colnames(smallRNAseqSamples) <- c("Sample", "Age", "Sex")
smallRNAseqSamples$Sample <- paste0("sample_", 1:nrow(smallRNAseqSamples))
write.csv(smallRNAseqSamples, "./data-raw/smallRNAseq_heart_samples.csv",
  row.names = FALSE
)
smallRNAseq_heart_samples <- smallRNAseqSamples
usethis::use_data(smallRNAseq_heart_samples, overwrite = TRUE)



# Proteomics data
protein10wks <- readxl::read_excel("./data-raw/Additional_file9.xlsx",
  sheet = 1
)
protein18wks <- readxl::read_excel("./data-raw/Additional_file10.xlsx",
  sheet = 1
)
colnames(protein10wks) <- c("gene", paste0("sample_", 1:3))
colnames(protein18wks) <- c("gene", paste0("sample_", 4:6))
protein <- full_join(protein10wks, protein18wks, by = "gene")
protein[is.na(protein)] <- 0
write.csv(protein, "./data-raw/protein_heart.csv", row.names = FALSE)
protein_heart <- protein
usethis::use_data(protein_heart, overwrite = TRUE)



# Sample matrix for proteomics data
proteinSamples <- data.frame(
  Sample = paste0("sample_", 1:6),
  Age = c(rep(10, 3), rep(18, 3))
)
write.csv(proteinSamples, "./data-raw/protein_heart_samples.csv",
  row.names = FALSE
)
protein_heart_samples <- proteinSamples
usethis::use_data(protein_heart_samples, overwrite = TRUE)


# Simulation of ATACseq peaks using some differential expressed genes
deResult <- read.csv("./data-raw/diff_expr_result_rnaseq.csv",
  header = TRUE, row.names = 1
)
# Retrieve up and down regulated genes
upGenes <- rownames(deResult)[deResult$logFC > 0 & deResult$padj < 0.05]
downGenes <- rownames(deResult)[deResult$logFC < 0 & deResult$padj < 0.05]
# Rank the genes
deResult$piVal <- log(deResult$padj, 10) * -1 * sign(deResult$logFC)
deResult <- deResult[order(deResult$piVal), ]
# Retrieve the top 100 up and down regulated genes
topUpGenes <- rownames(deResult)[1:200]
topDownGenes <- rownames(deResult)[(nrow(deResult) - 199):nrow(deResult)]

# Annotate TSS for all differentially expressed genes
library("biomaRt")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
GList <- getBM(
  filters = "hgnc_symbol", attributes =
    c(
      "hgnc_symbol", "start_position",
      "end_position", "strand", "chromosome_name"
    ),
  values = c(upGenes, downGenes), mart = mart
)

# Filter out genes that are not on the chromosomes 1-22, X, Y
GList <- GList[GList$chromosome_name %in% c(1:22, "X", "Y"), ]
# Add "chr" to the chromosome names
GList$chromosome_name <- paste0("chr", GList$chromosome_name)

# Simulating peaks as a general background for all DE genes
# Using the ATACseq data from GEO: GSE227923
# Randomly select 2000 peaks from the data
ATACseq <- read.table("./data-raw/GSM7111005_WT_ATAC_k_CM_rep1.narrowPeak",
  header = FALSE
)
# Retain only peaks from chr1-22, X, Y
ATACseq <- ATACseq[ATACseq$V1 %in% paste0("chr", c(1:22, "X", "Y")), ]
ATACseq <- ATACseq[, 1:3]
colnames(ATACseq) <- c("chr", "start", "end")

set.seed(91711)
# Randomly select 500 peaks
randomPeaks <- ATACseq[sample(1:nrow(ATACseq), 500), ]
set.seed(NULL)

# Simulate peaks TSS-based peaks for the top 200 up and down regulated genes
numPeaks <- 2000
peakWidth <- 500
simulatedPeaksUp <- data.frame(
  chr = character(),
  start = numeric(),
  end = numeric(),
  gene = character()
)
simulatedPeaksDown <- data.frame(
  chr = character(),
  start = numeric(),
  end = numeric(),
  gene = character()
)

topUpGenes <- topUpGenes[topUpGenes %in% GList$hgnc_symbol]
topDownGenes <- topDownGenes[topDownGenes %in% GList$hgnc_symbol]

set.seed(91711)
# Simulate peaks for up regulated genes
for (i in 1:numPeaks) {
  gene <- sample(topUpGenes, 1)
  tss <- GList[GList$hgnc_symbol == gene, "start_position"]
  strand <- GList[GList$hgnc_symbol == gene, "strand"]
  # TSS-proximal region defined to be 1000bp upstream and 3000bp downstream
  # of the TSS
  if (strand == 1) {
    promotorStart <- tss - 1500
    promotorEnd <- tss + 3000
  } else {
    promotorStart <- tss - 3000
    promotorEnd <- tss + 1500
  }
  # Randomly select a position within the promotor region
  peakStart <- sample(promotorStart:promotorEnd, 1)
  peakEnd <- peakStart + peakWidth
  # Append the simulated peak to the data frame
  simulatedPeaksUp <- rbind(
    simulatedPeaksUp,
    data.frame(
      chr = GList[
        GList$hgnc_symbol == gene,
        "chromosome_name"
      ],
      start = peakStart,
      end = peakEnd,
      gene = gene
    )
  )
}

# Simulate peaks for down regulated genes
for (i in 1:numPeaks) {
  gene <- sample(topDownGenes, 1)
  tss <- GList[GList$hgnc_symbol == gene, "start_position"]
  strand <- GList[GList$hgnc_symbol == gene, "strand"]
  # TSS-proximal region defined to be 1000bp upstream and 3000bp downstream
  # of the TSS
  if (strand == 1) {
    promotorStart <- tss - 1500
    promotorEnd <- tss + 3000
  } else {
    promotorStart <- tss - 3000
    promotorEnd <- tss + 1500
  }
  # Randomly select a position within the promotor region
  peakStart <- sample(promotorStart:promotorEnd, 1)
  # randomize the peak width with a size factor between 0.7 and 1.5
  peakEnd <- peakStart + peakWidth
  # Append the simulated peak to the data frame
  simulatedPeaksDown <- rbind(
    simulatedPeaksDown,
    data.frame(
      chr = GList[
        GList$hgnc_symbol == gene,
        "chromosome_name"
      ],
      start = peakStart,
      end = peakEnd,
      gene = gene
    )
  )
}
set.seed(NULL)


# Validate the size of the simulated peaks
simulatedPeaksUp$peakWidth <- simulatedPeaksUp$end - simulatedPeaksUp$start
summary(simulatedPeaksUp$peakWidth)
simulatedPeaksDown$peakWidth <- simulatedPeaksDown$end -
  simulatedPeaksDown$start
summary(simulatedPeaksDown$peakWidth)

# Save the simulated peaks
write.csv(simulatedPeaksUp, "./data-raw/simulatedPeaksUp.csv",
  row.names = FALSE
)
write.csv(simulatedPeaksDown, "./data-raw/simulatedPeaksDown.csv",
  row.names = FALSE
)

# Exchange peaks between up and down regulated genes with 100 randomly selected peaks
set.seed(91711)
selUp <- sample(1:nrow(simulatedPeaksUp), 100)
nselUp <- (1:nrow(simulatedPeaksUp))[-selUp]
selDown <- sample(1:nrow(simulatedPeaksDown), 100)
nselDown <- (1:nrow(simulatedPeaksDown))[-selDown]
peaksCond1 <- simulatedPeaksUp[selUp, ] %>% dplyr::select(chr, start, end)
peaksCond1 <- rbind(peaksCond1, simulatedPeaksDown[nselDown, ] %>% dplyr::select(chr, start, end))
peaksCond2 <- simulatedPeaksDown[selDown, ] %>% dplyr::select(chr, start, end)
peaksCond2 <- rbind(peaksCond2, simulatedPeaksUp[nselUp, ] %>% dplyr::select(chr, start, end))

# Add the random peaks to the simulated peaks, distribute them randomly
# to each condition
fraction <- runif(1, min = 0.3, max = 0.7)
numPeaksSel <- round(nrow(randomPeaks) * fraction)
selPeaks <- sample(1:nrow(randomPeaks), numPeaksSel)
nselPeaks <- (1:nrow(randomPeaks))[-selPeaks]
peaksCond1 <- rbind(peaksCond1, randomPeaks[selPeaks, ] %>%
  dplyr::select(chr, start, end))
peaksCond2 <- rbind(peaksCond2, randomPeaks[nselPeaks, ] %>%
  dplyr::select(chr, start, end))
set.seed(NULL)

# Save the simulated peaks
write.csv(peaksCond1, "./data-raw/peaksCond1.csv", row.names = FALSE)
write.csv(peaksCond2, "./data-raw/peaksCond2.csv", row.names = FALSE)

# Save to BED format
write.table(peaksCond1, "./data-raw/peaksCond1.bed",
  sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = FALSE
)
write.table(peaksCond2, "./data-raw/peaksCond2.bed",
  sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = FALSE
)



# Simply obtain the promoter regions of differentially expressed genes
# for motif enrichment analysis
peak1 <- data.frame(
  chr = character(),
  start = numeric(),
  end = numeric(),
  gene = character()
)
for (gene in topDownGenes) {
  tss <- GList[GList$hgnc_symbol == gene, "start_position"]
  strand <- GList[GList$hgnc_symbol == gene, "strand"]
  if (strand == 1) {
    str1 <- tss - 1000
    str2 <- tss - 500
    str3 <- tss
  } else {
    str1 <- tss + 500
    str2 <- tss
    str3 <- tss - 500
  }
  peak1 <- rbind(
    peak1,
    data.frame(
      chr = rep(GList[
        GList$hgnc_symbol == gene,
        "chromosome_name"
      ], 3),
      start = c(str1, str2, str3),
      end = c(str1 + 498, str2 + 498, str3 + 498),
      gene = rep(gene, 3)
    )
  )
}

peak2 <- data.frame(
  chr = character(),
  start = numeric(),
  end = numeric(),
  gene = character()
)
for (gene in topUpGenes) {
  tss <- GList[GList$hgnc_symbol == gene, "start_position"]
  strand <- GList[GList$hgnc_symbol == gene, "strand"]
  if (strand == 1) {
    str1 <- tss - 1000
    str2 <- tss - 500
    str3 <- tss
  } else {
    str1 <- tss + 500
    str2 <- tss
    str3 <- tss - 500
  }
  peak2 <- rbind(
    peak2,
    data.frame(
      chr = rep(GList[
        GList$hgnc_symbol == gene,
        "chromosome_name"
      ], 3),
      start = c(str1, str2, str3),
      end = c(str1 + 498, str2 + 498, str3 + 498),
      gene = rep(gene, 3)
    )
  )
}

peak1 <- peak1[, 1:3]
peak2 <- peak2[, 1:3]

# Add random peaks to the promoter regions
peak1 <- rbind(peak1, randomPeaks[selPeaks, ] %>%
  dplyr::select(chr, start, end))
peak2 <- rbind(peak2, randomPeaks[nselPeaks, ] %>%
  dplyr::select(chr, start, end))


write.table(peak1, "./data-raw/peak1.bed",
  sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = FALSE
)
write.table(peak2, "./data-raw/peak2.bed",
  sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = FALSE
)


# External interaction data
upregmiR2gene <- read.csv("./data-raw/upreg_mrna_gene2mir.csv", header = TRUE)
downregmiR2gene <- read.csv("./data-raw/downreg_mrna_gene2mir.csv", header = TRUE)
upregTF2gene <- read.csv("./data-raw/upreg_mrna_gene2tf.csv", header = TRUE)
downregTF2gene <- read.csv("./data-raw/downreg_mrna_gene2tf.csv", header = TRUE)

upregmiR2gene <- list(
  regulator = upregmiR2gene$regulator,
  target = upregmiR2gene$target
)
downregmiR2gene <- list(
  regulator = downregmiR2gene$regulator,
  target = downregmiR2gene$target
)
upregTF2gene <- list(
  regulator = upregTF2gene$regulator,
  target = upregTF2gene$target
)
downregTF2gene <- list(
  regulator = downregTF2gene$regulator,
  target = downregTF2gene$target
)

usethis::use_data(upregmiR2gene, downregmiR2gene, upregTF2gene,
  downregTF2gene,
  overwrite = TRUE
)


# Serialize JASPAR position weight matrices
library("JASPAR2022")
jasparVertebratePWM <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022,
  opts = list(
    matrixtype = "PWM",
    tax_group = "vertebrates"
  )
)


usethis::use_data(jasparVertebratePWM, overwrite = TRUE)


# Down sample RNAseq and small RNAseq data
# RNAseq
# Randomly select 5000 genes
set.seed(91711)
selGenes <- sample(1:nrow(RNAseq_heart), 5000)
rna <- RNAseq_heart[selGenes, ]
RNAseq_heart <- rna
usethis::use_data(RNAseq_heart, overwrite = TRUE)

# Small RNAseq
# Randomly select 5000 genes
# There are too many circRNAs in the small RNAseq data, so we will
# randomly select 1000 circRNAs and 5000 other genes
circIndex <- which(rownames(smallRNAseq_heart) %in% SNCANNOLIST_HSAPIENS$circRNA)
otherIndex <- which(!rownames(smallRNAseq_heart) %in% SNCANNOLIST_HSAPIENS$circRNA)
selCirc <- sample(circIndex, 100)
selOther <- sample(otherIndex, 0.7 * length(otherIndex))
selGenes <- c(selCirc, selOther)
srna <- smallRNAseq_heart[selGenes, ]
smallRNAseq_heart <- srna
usethis::use_data(smallRNAseq_heart, overwrite = TRUE)



# miR2Gene interaction data
miR2Gene <- read.csv("./data-raw/gene2mir.csv", header = TRUE)
set.seed(91711)
mir2Gene <- miR2Gene[sample(1:nrow(miR2Gene), 2000), ]
