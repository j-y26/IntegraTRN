# Purpose: Testing DETag
# Author: Jielin Yang
# Date: 2023-11-13
# Version: 1.0
# Bugs and Issues: None


# validateDETagSlots

test_that("validateDETagSlots throw errors if not conform to definition", {
  expect_error(validateDETagSlots(matrix(1:10), "DESeq2", NULL))
  expect_error(validateDETagSlots(data.frame(), "DESeq2", "something"))
  expect_error(validateDETagSlots(data.frame(), "something wrong", NULL))
})


# Constructor

# Create example data
# Define usable objects
deResult <- data.frame(
  log2FoldChange = c(1, 2, 3, -4, -5),
  pvalue = c(0.0001, 0.0002, 0.0003, 0.0004, 0.0005),
  padj = c(0.0004, 0.0005, 0.0006, 0.0007, 0.0008)
)
rownames(deResult) <- c("gene1", "gene2", "gene3", "gene4", "gene5")

deResult2 <- data.frame(
  logFC = c(1, 2, 3, -4, -5),
  PValue = c(0.0001, 0.0002, 0.0003, 0.0004, 0.0005),
  FDR = c(0.0004, 0.0005, 0.0006, 0.0007, 0.0008)
)
rownames(deResult2) <- c("gene1", "gene2", "gene3", "gene4", "gene5")

normCounts <- matrix(1:25, nrow = 5, ncol = 5)

peaks <- data.frame(
  chr = c("chr1", "chr2", "chr3", "chr4", "chr5"),
  start = c(1, 2, 3, 4, 5),
  end = c(2, 3, 4, 5, 6),
  Condition = c("-", "-", "+", "+", "+")
)

test_that("works for DESeq2 results", {
  detag <- DETag(deResult, "DESeq2", normCounts)
  expect_equal(detag@method, "DESeq2")
  expect_equal(detag@normalizedCounts, normCounts)
  expect_equal(detag@DEResult, deResult)
})

test_that("works for edgeR results", {
  detag <- DETag(deResult2, "edgeR", normCounts)
  expect_equal(detag@method, "edgeR")
  expect_equal(detag@normalizedCounts, normCounts)
  expect_equal(detag@DEResult, deResult2)
})

test_that("throws error when normalized count is null for count data", {
  expect_error(DETag(deResult, "DESeq2"))
})

test_that("works for ATAC data", {
  detag <- DETag(peaks, "GRangesATAC")
  expect_equal(detag@method, "GRangesATAC")
  expect_equal(length(detag@normalizedCounts), 0)
  expect_equal(detag@DEResult, peaks)
})

test_that("throws error when method is invalid", {
  expect_error(DETag(deResult, "something wrong", normCounts))
})


# exportDE

test_that("exportDE exports defined columns when original = FALSE", {
  # DESeq2
  detag <- DETag(deResult, "DESeq2", normCounts)
  deres <- exportDE(detag, original = FALSE)
  expect_equal(colnames(deres), COUNT_DEFIELDS)
  expect_true(all(rownames(deres) %in% rownames(deResult)))
  # edgeR
  detag <- DETag(deResult2, "edgeR", normCounts)
  deres <- exportDE(detag, original = FALSE)
  expect_equal(colnames(deres), COUNT_DEFIELDS)
  expect_true(all(rownames(deres) %in% rownames(deResult2)))
  # ATAC
  detag <- DETag(peaks, "GRangesATAC")
  deres <- exportDE(detag, original = FALSE)
  expect_equal(colnames(deres), colnames(peaks))
})


test_that("exportDE exports original columns when original = TRUE", {
  # DESeq2
  detag <- DETag(deResult, "DESeq2", normCounts)
  deres <- exportDE(detag, original = TRUE)
  expect_equal(colnames(deres), DESEQ2_FIELDS)
  expect_true(all(rownames(deres) %in% rownames(deResult)))
  # edgeR
  detag <- DETag(deResult2, "edgeR", normCounts)
  deres <- exportDE(detag, original = TRUE)
  expect_equal(colnames(deres), EDGER_FIELDS)
  expect_true(all(rownames(deres) %in% rownames(deResult2)))
  # ATAC
  detag <- DETag(peaks, "GRangesATAC")
  deres <- exportDE(detag, original = TRUE)
  expect_equal(colnames(deres), colnames(peaks))
})


test_that("throws error if original not boolean", {
  detag <- DETag(deResult, "DESeq2", normCounts)
  expect_error(exportDE(detag, original = "something wrong"))
})


# exportNormalizedCounts

test_that("exportNormalizedCounts exports normalized counts", {
  detag <- DETag(deResult, "DESeq2", normCounts)
  expect_equal(colnames(exportNormalizedCounts(detag)), colnames(normCounts))
  expect_equal(rownames(exportNormalizedCounts(detag)), rownames(normCounts))
})

test_that("throws error when trying to export it of ATAC", {
  detag <- DETag(peaks, "GRangesATAC")
  expect_error(exportNormalizedCounts(detag))
})

# [END]
