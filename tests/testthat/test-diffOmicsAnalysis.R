# Purpose: Testing differential omics analysis
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


# filterGeneCounts
test_that("filterGeneCounts should not filter anything", {
  mol1 <- filterGeneCounts(mol, "RNAseq")
  expect_true(all(mol1@RNAseq == mol@RNAseq))
  expect_true(all(mol1@smallRNAseq == mol@smallRNAseq))
  mol1 <- filterGeneCounts(mol, "smallRNAseq")
  expect_true(all(mol1@RNAseq == mol@RNAseq))
  expect_true(all(mol1@smallRNAseq == mol@smallRNAseq))
  expect_error(filterGeneCounts(mol, "somethingwrong"))
})

test_that("filter something that do not exist", {
  expect_silent(filterGeneCounts(mol, "proteomics"))
  mol1 <- filterGeneCounts(mol, "proteomics")
  expect_true(all(mol1@RNAseq == mol@RNAseq))
})

test_that("filterout zero counts", {
  mol1 <- mol
  mol1@RNAseq[8:10, 1:9] <- 0
  mol1 <- filterGeneCounts(mol1, "RNAseq")
  expect_equal(nrow(mol1@RNAseq), 7)
  expect_equal(ncol(mol1@RNAseq), 10)
  expect_equal(nrow(mol1@smallRNAseq), 10)
  expect_equal(ncol(mol1@smallRNAseq), 10)
  mol1 <- filterGeneCounts(mol1, "smallRNAseq")
  expect_equal(nrow(mol1@RNAseq), 7)
  expect_equal(ncol(mol1@RNAseq), 10)
  expect_equal(nrow(mol1@smallRNAseq), 10)
  expect_equal(ncol(mol1@smallRNAseq), 10)
})


# diffOmics

test_that("test diffOmics runs and ignores non-exiting omics", {
  mol1 <- diffOmics(mol, program = "edgeR")
  expect_false(is.null(mol1$DERNAseq))
  expect_false(is.null(mol1$DEsmallRNAseq))
  expect_true(is.null(mol1$DEproteomics))
  # Test type creation
  expect_s4_class(mol1$DERNAseq, "DETag")
  expect_s4_class(mol1$DEsmallRNAseq, "DETag")
})


# [END]
