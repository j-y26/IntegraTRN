# Purpose: Testing processSmallRNA
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

# Sample annotation
sampleAnno <- list(
  miRNA = paste0("transcript", 1:3),
  piRNA = paste0("transcript", 4:6),
  tRNA = paste0("transcript", 7:10)
)
mol$annoSncRNA <- sampleAnno


# sncAnnoCoverage

test_that("invalid anno param", {
  expect_error(sncAnnoCoverage(allTranscripts, mol, anno = "invalid"))
})

test_that("full annotation", {
  passed <- sncAnnoCoverage(allTranscripts, mol, anno = "userAnno")
  expect_true(passed)
})

test_that("partial annotation", {
  transcripts <- paste0("transcript", 1:20)
  passed <- sncAnnoCoverage(transcripts, mol, anno = "userAnno")
  expect_false(passed)
})

# checkSmallAnnoCoverage

test_that("invalid anno param", {
  expect_error(checkSmallAnnoCoverage(mol, anno = "invalid"))
})

test_that("full annotation", {
  expect_silent(checkSmallAnnoCoverage(mol, anno = "userAnno"))
})


# extractTranscriptFromAnno

# Create an annoDF
annoDF <- data.frame(
  transcript = paste0("transcript", 1:10),
  type = c(rep("miRNA", 3), rep("piRNA", 3), rep("tRNA", 4))
)

test_that("invalid category", {
  expect_warning(extractTranscriptFromAnno(annoDF, "invalid"))
})

test_that("category at start of DF", {
  expect_true(all(extractTranscriptFromAnno(annoDF, "miRNA") %in% paste0("transcript", 1:3)))
})

test_that("category at end of DF", {
  expect_true(all(extractTranscriptFromAnno(annoDF, "tRNA") %in% paste0("transcript", 7:10)))
})

test_that("category in middle of DF", {
  expect_true(all(extractTranscriptFromAnno(annoDF, "piRNA") %in% paste0("transcript", 4:6)))
})

# annotateSmallRNA

test_that("invalid anno param", {
  expect_error(annotateSmallRNA(mol, anno = "invalid"))
})

# make the annotation cover all types
annoDF <- rbind(
  annoDF,
  data.frame(
    transcript = paste0("transcript", 11:13),
    type = c("snoRNA", "snRNA", "circRNA")
  )
)

test_that("full annotation", {
  expect_silent(annotateSmallRNA(mol, anno = annoDF))
})


# [END]
