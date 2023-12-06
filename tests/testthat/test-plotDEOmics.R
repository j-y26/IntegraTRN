# Purpose: Testing plotting functions
# Author: Jielin Yang
# Date: 2023-12-06
# Version: 1.0
# Bugs and Issues: None


# Define basic data for testing
rnaseqCountMatrix <- matrix(sample(0:100, 1000, replace = TRUE),
  nrow = 100, ncol = 10
)
rownames(rnaseqCountMatrix) <- paste0("Gene", 1:100)
colnames(rnaseqCountMatrix) <- paste0("Sample", 1:10)
rnaseqGroupBy <- rep(c("Group1", "Group2"), each = 5)

smallRNAseqCountMatrix <- matrix(sample(0:100, 800, replace = TRUE),
  nrow = 100, ncol = 8
)
rownames(smallRNAseqCountMatrix) <- paste0("Gene", 1:100)
colnames(smallRNAseqCountMatrix) <- paste0("Sample", 1:8)
smallRNAseqGroupBy <- rep(c("Group1", "Group2"), each = 4)

proteomicsCountMatrix <- matrix(sample(0:100, 600, replace = TRUE),
  nrow = 100, ncol = 6
)
rownames(proteomicsCountMatrix) <- paste0("Gene", 1:100)
colnames(proteomicsCountMatrix) <- paste0("Sample", 1:6)
proteomicsGroupBy <- rep(c("Group1", "Group2"), each = 3)

# Randomly select some peaks from the package data
peak1 <- system.file("extdata", "peak1.bed", package = "IntegraTRN")
peak2 <- system.file("extdata", "peak2.bed", package = "IntegraTRN")
peak1 <- read.table(peak1, header = FALSE, sep = "\t")
peak2 <- read.table(peak2, header = FALSE, sep = "\t")
peak1 <- peak1[sample(1:nrow(peak1), 100), ]
peak2 <- peak2[sample(1:nrow(peak2), 100), ]

# Integration tests need to be run for testing the plotting functions

test_that("test plotVolcano generates correct ggplot", {
  # Create a MOList
  objMOList <- MOList(
    RNAseq = rnaseqCountMatrix, RNAGroupBy = rnaseqGroupBy
  )

  # No differential expression analysis, should render error on plotting
  expect_error(plotVolcano(objMOList)) # default
  expect_error(plotVolcano(objMOList, "RNAseq")) # specify data type
  expect_error(plotVolcano(objMOList, "smallRNAseq")) # specify group
  expect_error(plotVolcano(objMOList, "proteomics")) # specify group
  expect_error(plotVolcano(objMOList, "something")) # wrong group

  # Add DE analysis
  objMOList <- diffOmics(objMOList, program = "edgeR")

  # Plotting
  expect_silent(plotVolcano(objMOList)) # default
  rnaPlt <- plotVolcano(objMOList, "RNAseq") # specify data type
  expect_s3_class(rnaPlt, "gg")
  expect_true(is.null(rnaPlt$labels$title))
  expect_error(plotVolcano(objMOList, "smallRNAseq")) # no data
  expect_error(plotVolcano(objMOList, "proteomics")) # no data
  rnaPlt1 <- plotVolcano(objMOList, "RNAseq", title = "someTitle")
  expect_s3_class(rnaPlt1, "gg")
  expect_equal(rnaPlt1$labels$title, "someTitle")

  # Create a MOList with all data
  objMOList <- MOList(
    RNAseq = rnaseqCountMatrix, RNAGroupBy = rnaseqGroupBy,
    smallRNAseq = smallRNAseqCountMatrix, smallRNAGroupBy = smallRNAseqGroupBy,
    proteomics = proteomicsCountMatrix, proteomicsGroupBy = proteomicsGroupBy,
    ATACpeak1 = peak1, ATACpeak2 = peak2
  )
  objMOList <- diffOmics(objMOList, program = "edgeR")

  # Plotting
  expect_silent(plotVolcano(objMOList)) # default
  rnaPlt <- plotVolcano(objMOList, "RNAseq") # specify data type
  expect_s3_class(rnaPlt, "gg")
  expect_true(is.null(rnaPlt$labels$title))
  expect_silent(plotVolcano(objMOList, "smallRNAseq")) # specify group
  expect_silent(plotVolcano(objMOList, "proteomics")) # specify group
  expect_error(plotVolcano(objMOList, "something")) # wrong group
  smallRnaPlt <- plotVolcano(objMOList, "smallRNAseq")
  expect_s3_class(smallRnaPlt, "gg")
  expect_true(is.null(smallRnaPlt$labels$title))
  proteomicsPlt <- plotVolcano(objMOList, "proteomics")
  expect_s3_class(proteomicsPlt, "gg")
  expect_true(is.null(proteomicsPlt$labels$title))
})

test_that("test plotVolcanoSmallRNA generates correct outputs", {
  # Create a MOList
  objMOList <- MOList(
    RNAseq = rnaseqCountMatrix, RNAGroupBy = rnaseqGroupBy,
    smallRNAseq = smallRNAseq_heart[1:100, 1:8],
    smallRNAGroupBy = smallRNAseqGroupBy
  )

  # No differential expression analysis, should render error on plotting
  expect_error(plotVolcanoSmallRNA(objMOList)) # default
  expect_error(plotVolcanoSmallRNA(objMOList, "smallRNAseq")) # specify data type
  expect_error(plotVolcanoSmallRNA(objMOList, "RNAseq")) # specify group
  expect_error(plotVolcanoSmallRNA(objMOList, "proteomics")) # specify group
  expect_error(plotVolcanoSmallRNA(objMOList, "something")) # wrong group

  # Add DE analysis
  objMOList <- diffOmics(objMOList, program = "edgeR")

  # No annotation, should render error on plotting
  expect_error(plotVolcanoSmallRNA(objMOList))

  # Add annotation
  objMOList <- annotateSmallRNA(objMOList)

  # Plotting
  expect_silent(plotVolcanoSmallRNA(objMOList)) # default
  smallRnaPlt <- plotVolcanoSmallRNA(objMOList) # default
  expect_s3_class(smallRnaPlt, "gg")
  expect_true(is.null(smallRnaPlt$labels$title))
  smallRnaPlt1 <- plotVolcanoSmallRNA(objMOList, title = "someTitle")
  expect_s3_class(smallRnaPlt1, "gg")
  expect_equal(smallRnaPlt1$labels$title, "someTitle")
})

test_that("test countPCA", {
  # Create a non-numeric matrix
  failMatrix <- matrix(sample(letters, 100, replace = TRUE),
    nrow = 10, ncol = 10
  )
  group <- rep(c("Group1", "Group2"), each = 5)
  expect_error(countPCA(failMatrix, group), "numeric matrix")

  # one-group groupBy
  expect_error(countPCA(rnaseqCountMatrix, rep("g1", 10)), "two groups")

  # wrong groupBy length
  group <- rep(c("Group1", "Group2"), each = 4)
  expect_error(countPCA(rnaseqCountMatrix, group), "equal")

  # batch var the same
  batch <- rep("Batch1", 10)
  batch <- as.factor(batch)
  rnaseqGroupBy <- as.factor(rnaseqGroupBy)
  expect_warning(countPCA(rnaseqCountMatrix, rnaseqGroupBy, batch))
  batch <- rep("Batch1", 5)
  batch <- as.factor(batch)
  expect_warning(countPCA(rnaseqCountMatrix, rnaseqGroupBy, batch))
  batch <- c("Batch1", "Batch2")
  batch <- as.factor(batch)
  expect_error(countPCA(rnaseqCountMatrix, rnaseqGroupBy, batch), "equal")

  # Valid input
  batch <- rep(c("Batch1", "Batch2", "Batch3", "Batch4", "Batch5"), 2)
  batch <- as.factor(batch)
  rnaseqGroupBy <- as.factor(rnaseqGroupBy)
  expect_message(countPCA(rnaseqCountMatrix, rnaseqGroupBy, batch))
})

test_that("plotSmallRNAPCAs works correctlys", {
  # Create a MOList
  objMOList <- MOList(
    RNAseq = rnaseqCountMatrix, RNAGroupBy = as.factor(rnaseqGroupBy),
    smallRNAseq = smallRNAseq_heart[1:100, 1:8],
    smallRNAGroupBy = as.factor(smallRNAseqGroupBy)
  )

  # No annotation, should render error on plotting
  expect_error(plotSmallRNAPCAs(objMOList))

  # Add annotation
  objMOList <- annotateSmallRNA(objMOList)

  # Plotting
  expect_message(plotSmallRNAPCAs(objMOList)) # default
  smallRnaPlts <- plotSmallRNAPCAs(objMOList) # default
  expect_true(is.list(smallRnaPlts))
  expect_s3_class(smallRnaPlts[[1]], "gg")
})


test_that("plotting functions for atacseq data", {
  # Create a MOList
  objMOList <- MOList(
    RNAseq = rnaseqCountMatrix, RNAGroupBy = rnaseqGroupBy,
    ATACpeak1 = peak1, ATACpeak2 = peak2
  )

  # No annotation, should render error on plotting
  expect_error(plotATACCoverage(objMOList))
  expect_error(plotATACAnno(objMOList))
  expect_error(plotATACMotifHeatmap(objMOList))

  # differential analysis, but not for peaks
  objMOList <- diffOmics(objMOList, program = "edgeR")
  expect_silent(plotATACCoverage(objMOList)) # can be plotted
  expect_error(plotATACAnno(objMOList))
  expect_error(plotATACMotifHeatmap(objMOList))

  # motif annotation and enrichment
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  annodb <- "org.Hs.eg.db"
  bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  jasparVertebratePWM <- jasparVertebratePWM[1:10]
  objMOList <- annotateATACPeaksMotif(
    objMOList,
    tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = annodb,
    bsgenome = bsgenome, pwmL = jasparVertebratePWM
  )
  expect_silent(plotATACCoverage(objMOList))
  expect_silent(plotATACAnno(objMOList))
})
