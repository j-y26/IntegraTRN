# Purpose: Testing ATAC Peak Processing
# Author: Jielin Yang
# Date: 2023-11-13
# Version: 1.0
# Bugs and Issues: None

# mergePeaks

# Create an example BED dataframe
test_that("no overlap", {
  peaks <- data.frame(
    chr = c("chr1", "chr1", "chr2"),
    start = c(1, 5, 9),
    end = c(2, 6, 10)
  )
  peakGR <- mergePeaks(peaks)
  expect_equal(length(peakGR), 3)
})

test_that("overlaps", {
  peaks <- data.frame(
    chr = c("chr1", "chr1", "chr2"),
    start = c(1, 2, 9),
    end = c(2, 6, 10)
  )
  peakGR <- mergePeaks(peaks)
  expect_equal(length(peakGR), 2)
})


# Create an example MOList
rnamatrix <- matrix(1:100, nrow = 10, ncol = 10)
rownames(rnamatrix) <- paste0("gene", 1:10)
colnames(rnamatrix) <- paste0("sample", 1:10)
rnaGroupBy <- c(
  "group1", "group1", "group1", "group1", "group1",
  "group2", "group2", "group2", "group2", "group2"
)
smallrnamatrix <- matrix(1:100, nrow = 10, ncol = 10)
allTranscripts <- paste0("transcript", 1:10)
rownames(smallrnamatrix) <- allTranscripts
colnames(smallrnamatrix) <- paste0("sample", 1:10)
smallrnaGroupBy <- c(
  "group1", "group1", "group1", "group1", "group1",
  "group2", "group2", "group2", "group2", "group2"
)
mol <- MOList(
  RNAseq = rnamatrix, RNAGroupBy = rnaGroupBy,
  smallRNAseq = smallrnamatrix, smallRNAGroupBy = smallrnaGroupBy
)

test_that("no DEATAC", {
  bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  data("jasparVertebratePWM", package = "IntegraTRN")
  expect_error(enrichMotifs(mol, bsgenome, jasparVertebratePWM))
})

test_that("invalid pwmlist", {
  bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  expect_error(enrichMotifs(mol, bsgenome, "something invalid"))
})

# Integration test on simulated data exported by package IntegraTRN
# Very hard to design light-weight unit tests
data("expMOList", package = "IntegraTRN") # contains some ATAC peaks
data("jasparVertebratePWM", package = "IntegraTRN")
test_that("annotate and enrich for motifs", {
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  annodb <- "org.Hs.eg.db"
  # annotate peaks
  expMOList <- annotatePeaks(expMOList, c(-3000, 3000), txdb, annodb)
  expect_false(is.null(expMOList$DEATAC))
  expect_true(inherits(expMOList$DEATAC, "PEAKTag"))
  expect_equal(class(expMOList$DEATAC@TxDB)[1], "TxDb")
  expect_equal(expMOList$DEATAC@annoDB, "org.Hs.eg.db")
  expect_equal(class(expMOList$DEATAC@annotatedPeaks)[1], "csAnno")

  bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  # motif enrichment
  expMOList <- enrichMotifs(expMOList, bsgenome, jasparVertebratePWM)
  expect_false(is.null(expMOList$DEATAC$motifEnrichment))
  expect_equal(class(expMOList$DEATAC$motifEnrichment)[1], "SummarizedExperiment")
})


# [END]
