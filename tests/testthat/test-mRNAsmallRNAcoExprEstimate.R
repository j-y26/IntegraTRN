# Purpose: Testing xo-expression estimation between mRNA and smallRNA
# Author: Jielin Yang
# Date: 2023-11-13
# Version: 1.0
# Bugs and Issues: None

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

# validateSampleDF

sampleDFRNA <- data.frame(
  age = 1:10,
  group = rnaGroupBy
)
sampleDFSmallRNA <- data.frame(
  age = 1:10,
  group = smallrnaGroupBy
)

test_that("df matches with molist", {
  rownames(sampleDFRNA) <- colnames(rnamatrix)
  rownames(sampleDFSmallRNA) <- colnames(smallrnamatrix)
  expect_silent(validateSampleDFs(mol, sampleDFRNA, sampleDFSmallRNA))
})

test_that("name mismatch", {
  expect_error(validateSampleDFs(mol, sampleDFRNA, sampleDFSmallRNA))
})

test_that("mismatch dimension", {
  sampleDFSmallRNA <- data.frame(
    age = 1:11,
    group = c(smallrnaGroupBy, "group2")
  )
  rownames(sampleDFSmallRNA) <- c(colnames(smallrnamatrix), "gene")
  expect_error(validateSampleDFs(mol, sampleDFRNA, sampleDFSmallRNA))
})

test_that("mismatch dimension another", {
  sampleDFRNA <- data.frame(
    age = 1:11,
    group = c(rnaGroupBy, "group2")
  )
  rownames(sampleDFRNA) <- c(colnames(rnamatrix), "gene")
  expect_error(validateSampleDFs(mol, sampleDFRNA, sampleDFSmallRNA))
})


# matchContinuous

# Create an example sampleDF with continuous grouping variable

sampleDF <- data.frame(
  age = rep(1:10, 2),
  groupBy = sample(1:5, 20, replace = TRUE)
)
rownames(sampleDF) <- paste0("sample", 1:20)
rownames(sampleDF)[13:20] <- paste0(rownames(sampleDF)[13:20], SRNA_SUFFIX)

test_that("matchContinuous works", {
  expect_silent(matchContinuous(sampleDF))
  results <- matchContinuous(sampleDF)
  # Test for result formatting
  expect_true(is.data.frame(results))
  expect_true("RNA" %in% colnames(results))
  expect_true("smallRNA" %in% colnames(results))
  expect_equal(ncol(results), 2)
})

# matchInGroup

# Create an example sampleDF with no grouping variable

sampleDF1 <- sampleDF[, 1, drop = FALSE]
sampleDF1$seq <- c(rep(0, 12), rep(1, 8))

test_that("matchInGroup works", {
  expect_silent(matchInGroup(sampleDF1))
  results <- matchInGroup(sampleDF1)
  # Test for result formatting
  expect_true(is.data.frame(results))
  expect_true("RNA" %in% colnames(results))
  expect_true("smallRNA" %in% colnames(results))
  expect_equal(ncol(results), 2)
})


# matchBinary

# Create an example sampleDF with binary grouping variable

sampleDF2 <- sampleDF
sampleDF2$groupBy <- sampleDF2$groupBy > 3

test_that("matchBinary works", {
  expect_silent(matchBinary(sampleDF2))
  results <- matchBinary(sampleDF2)
  # Test for result formatting
  expect_true(is.data.frame(results))
  expect_true("RNA" %in% colnames(results))
  expect_true("smallRNA" %in% colnames(results))
  expect_equal(ncol(results), 2)
})


# nnRNAMatch, which is upstream of the above functions

sampleDFRNA <- data.frame(
  age = 1:10,
  groupBy = rnaGroupBy
)
rownames(sampleDFRNA) <- colnames(rnamatrix)
sampleDFSmallRNA <- data.frame(
  age = 1:10,
  groupBy = smallrnaGroupBy
)
rownames(sampleDFSmallRNA) <- colnames(smallrnamatrix)

test_that("nnRNAMatch works with same ordered match", {
  result <- nnRNAMatch(mol, sampleDFRNA, sampleDFSmallRNA)
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  expect_equal(length(result$indexRNAseq), 10)
  expect_equal(length(result$indexSmallRNAseq), 10)
  expect_true(all(result$indexRNAseq == 1:10))
  expect_true(all(result$indexSmallRNAseq == 1:10))
})

test_that("nnRNA works with different ordered match", {
  sampleDFRNA <- sampleDFRNA[10:1, ]
  result <- nnRNAMatch(mol, sampleDFRNA, sampleDFSmallRNA)
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  expect_equal(length(result$indexRNAseq), 10)
  expect_equal(length(result$indexSmallRNAseq), 10)
})


# matchSamplesRNAsmallRNA, even more upstream

test_that("matchSamplesRNAsmallRNA works with no additional data", {
  mol1 <- matchSamplesRNAsmallRNA(mol)
  expect_false(is.null(mol1$matchingRNAsmallRNA))
  expect_equal(length(mol1$matchingRNAsmallRNA), 2)
  expect_equal(length(mol1$matchingRNAsmallRNA$indexRNAseq), 10)
  expect_equal(length(mol1$matchingRNAsmallRNA$indexSmallRNAseq), 10)
  # since the matches are random
})

test_that("matchSamplesRNAsmallRNA works with additional data", {
  mol1 <- matchSamplesRNAsmallRNA(mol, sampleDFRNA, sampleDFSmallRNA)
  expect_false(is.null(mol1$matchingRNAsmallRNA))
  expect_equal(length(mol1$matchingRNAsmallRNA), 2)
  expect_equal(length(mol1$matchingRNAsmallRNA$indexRNAseq), 10)
  expect_equal(length(mol1$matchingRNAsmallRNA$indexSmallRNAseq), 10)
  expect_true(all(mol1$matchingRNAsmallRNA$indexRNAseq == 1:10))
  expect_true(all(mol1$matchingRNAsmallRNA$indexSmallRNAseq == 1:10))
})


# exportMatchResult

test_that("export correct sample names", {
  mol1 <- matchSamplesRNAsmallRNA(mol, sampleDFRNA, sampleDFSmallRNA)
  result <- exportMatchResult(mol1)
  expect_true(is.data.frame(result))
  expect_equal(ncol(result), 2)
  expect_equal(nrow(result), 10)
  expect_true(all(result$RNAseq %in% colnames(rnamatrix)))
  expect_true(all(result$SmallRNAseq %in% colnames(smallrnamatrix)))
})

# [END]
