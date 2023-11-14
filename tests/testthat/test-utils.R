# Purpose: Testing untilities for the package
# Author: Jielin Yang
# Date: 2023-11-13
# Version: 1.0
# Bugs and Issues: None

# getSampleNames

test_that("matrix has no column names", {
  countMatrix <- matrix(1:9, nrow = 3, ncol = 3)
  expect_equal(getSampleNames(countMatrix), paste0("sample_", seq_len(ncol(countMatrix))))
})

test_that("matrix has column names", {
  countMatrix <- matrix(1:9, nrow = 3, ncol = 3)
  colnames(countMatrix) <- c("sample_1", "sample_2", "sample_3")
  expect_equal(getSampleNames(countMatrix), colnames(countMatrix))
})

test_that("different sample names than default", {
  countMatrix <- matrix(1:9, nrow = 3, ncol = 3)
  colnames(countMatrix) <- c("some", "random", "names")
  expect_equal(getSampleNames(countMatrix), colnames(countMatrix))
})

# matchVecToMatrix

test_that("vector and matrix have same length", {
  vec <- c("a", "b", "c")
  mat <- matrix(1:9, nrow = 3, ncol = 3)
  expect_true(matchVecToMatrix(vec, mat))
})

test_that("vector and matrix have different length", {
  vec <- c("a", "b", "c")
  mat <- matrix(1:6, nrow = 3, ncol = 2)
  expect_false(matchVecToMatrix(vec, mat))
})

# DESeqDesign

test_that("batch is NULL", {
  groupBy <- c("a", "b", "c")
  colDF <- data.frame(groupBy)
  colnames(colDF) <- "group"
  expect_equal(DESeqDesign(groupBy)$colData, colDF)
  expect_equal(class(DESeqDesign(groupBy)$design), "formula")
})

test_that("batch is not NULL", {
  groupBy <- c("a", "b", "c")
  batch <- c("x", "y", "z")
  colDF <- data.frame(groupBy, batch)
  colnames(colDF) <- c("group", "batch")
  expect_equal(DESeqDesign(groupBy, batch)$colData, colDF)
  expect_equal(class(DESeqDesign(groupBy, batch)$design), "formula")
})

# matchResultToDF

test_that("match result to data frame", {
  # Generate a example data frame for matching
  df <- data.frame(
    age = c(12, 14, 10),
    sex = c("M", "F", "M"),
    treat = c("A", "B", "A"),
    seq = c(1, 0, 0)
  )
  # Name the samples
  rownames(df) <- c("sample_1",
                    paste0("sample_2", SRNA_SUFFIX),
                    paste0("sample_3", SRNA_SUFFIX))
  # Generate a example result
  matchResult <- MatchIt::matchit(
    seq ~ age + sex + treat,
    data = df,
    method = "optimal",
    distance = "mahalanobis",
    ratio = 1
  )
  # Convert to data frame
  matchDF <- matchResultToDF(matchResult)
  # Check the result
  expect_equal(colnames(matchDF), c("smallRNA", "RNA"))
  expect_equal(nrow(matchDF), 1)
  expect_false(any(grepl(SRNA_SUFFIX, matchDF$smallRNA)))
  expect_false(any(grepl(SRNA_SUFFIX, matchDF$RNA)))
})


# labelSmallSizeGroup

test_that("RNA samples fewer", {
  df <- data.frame(
    age = c(12, 14, 10),
    sex = c("M", "F", "M"),
    treat = c("A", "B", "A")
  )
  # Name the samples
  rownames(df) <- c("sample_1",
                    paste0("sample_2", SRNA_SUFFIX),
                    paste0("sample_3", SRNA_SUFFIX))
  # Label with "seq" column
  df <- labelSmallSizeGroup(df, SRNA_SUFFIX, "seq")
  # Check the result
  expect_equal(df$seq, c(1, 0, 0))
})

test_that("smallrnaseq samples fewer", {
  df <- data.frame(
    age = c(12, 14, 10),
    sex = c("M", "F", "M"),
    treat = c("A", "B", "A")
  )
  # Name the samples
  rownames(df) <- c("sample_1",
                    "sample_2",
                    paste0("sample_3", SRNA_SUFFIX))
  # Label with "seq" column
  df <- labelSmallSizeGroup(df, SRNA_SUFFIX, "seq")
  # Check the result
  expect_equal(df$seq, c(0, 0, 1))
})

test_that("different colname", {
  df <- data.frame(
    age = c(12, 14, 10),
    sex = c("M", "F", "M"),
    treat = c("A", "B", "A")
  )
  # Name the samples
  rownames(df) <- c("sample_1",
                    paste0("sample_2", SRNA_SUFFIX),
                    paste0("sample_3", SRNA_SUFFIX))
  # Label with "seq" column
  df <- labelSmallSizeGroup(df, SRNA_SUFFIX, "something")
  # Check the result
  expect_equal(df$something, c(1, 0, 0))
})


# extractDirectionalGenes

test_that("no upregulated genes", {
  # fake DE results
  deRes <- data.frame(
      logFC = c(-1, -2, -3),
      padj = c(0.1, 0.2, 0.3)
  )
  rownames(deRes) <- c("gene_1", "gene_2", "gene_3")
  # extract directional genes
  dirGenes <- extractDirectionalGenes(deRes)
  # check the result
  expect_equal(dirGenes$up, character(0))
  expect_equal(dirGenes$down, c("gene_1", "gene_2", "gene_3"))
})

test_that("no downregulated genes", {
  # fake DE results
  deRes <- data.frame(
      logFC = c(1, 2, 3),
      padj = c(0.1, 0.2, 0.3)
  )
  rownames(deRes) <- c("gene_1", "gene_2", "gene_3")
  # extract directional genes
  dirGenes <- extractDirectionalGenes(deRes)
  # check the result
  expect_equal(dirGenes$up, c("gene_1", "gene_2", "gene_3"))
  expect_equal(dirGenes$down, character(0))
})

test_that("both upregulated and downregulated genes", {
  # fake DE results
  deRes <- data.frame(
      logFC = c(-1, 2, 3),
      padj = c(0.1, 0.2, 0.3)
  )
  rownames(deRes) <- c("gene_1", "gene_2", "gene_3")
  # extract directional genes
  dirGenes <- extractDirectionalGenes(deRes)
  # check the result
  expect_equal(dirGenes$up, c("gene_2", "gene_3"))
  expect_equal(dirGenes$down, c("gene_1"))
})

test_that("empty DE results", {
  # fake DE results
  deRes <- data.frame(
      logFC = numeric(0),
      padj = numeric(0)
  )
  rownames(deRes) <- character(0)
  # extract directional genes
  dirGenes <- extractDirectionalGenes(deRes)
  # check the result
  expect_equal(dirGenes$up, character(0))
  expect_equal(dirGenes$down, character(0))
})

# findGeneType

test_that("gene in the first type", {
  # fake annotation
  annotation <- list(
    type1 = c("gene_1", "gene_2"),
    type2 = c("gene_3", "gene_4")
  )
  # find gene type
  geneType <- findGeneType("gene_1", annotation)
  expect_true(geneType == "type1")
  # find gene type
  geneType <- findGeneType("gene_2", annotation)
  expect_true(geneType == "type1")
})

test_that("gene in the second type", {
  # fake annotation
  annotation <- list(
    type1 = c("gene_1", "gene_2"),
    type2 = c("gene_3", "gene_4")
  )
  # find gene type
  geneType <- findGeneType("gene_3", annotation)
  expect_true(geneType == "type2")
  # find gene type
  geneType <- findGeneType("gene_4", annotation)
  expect_true(geneType == "type2")
})

test_that("gene in no type", {
  # fake annotation
  annotation <- list(
    type1 = c("gene_1", "gene_2"),
    type2 = c("gene_3", "gene_4")
  )
  # find gene type
  geneType <- findGeneType("gene_5", annotation)
  expect_false(is.null(geneType["$gene_5"]))
})

test_that("empty annotation", {
  # fake annotation
  annotation <- list()
  # find gene type
  geneType <- findGeneType("gene_1", annotation)
  expect_false(is.null(geneType))
})

test_that("vectorized", {
  # fake annotation
  annotation <- list(
    type1 = c("gene_1", "gene_2"),
    type2 = c("gene_3", "gene_4")
  )
  # find gene type
  geneType <- findGeneType(c("gene_1", "gene_3"), annotation)
  expect_true(all(geneType == c("type1", "type2")))
})

# [END]