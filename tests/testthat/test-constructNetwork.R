# Purpose: Testing constructNetwork.R
# Author: Jielin Yang
# Date: 2023-12-07
# Version: 1.0
# Bugs and Issues: None


# Define example data
rnaCountMatrix <- matrix(1:100, nrow = 10, ncol = 10)
colnames(rnaCountMatrix) <- paste0("sample", 1:10)
rownames(rnaCountMatrix) <- paste0("gene", 1:10)
rnaSampleGroup <- c(rep("group1", 5), rep("group2", 5))

testMOList <- MOList(RNAseq = rnaCountMatrix, RNAGroupBy = rnaSampleGroup)

# External interaction data for testing
mir2genes <- list(
  regulator = c("mir1", "mir1", "mir2", "mir2", "mir3", "mir3"),
  target = c("gene1", "gene2", "gene3", "gene4", "gene5", "gene6")
)

tf2genes <- list(
  regulator = c("tf1", "tf1", "tf2", "tf2", "tf3", "tf3"),
  target = c("gene1", "gene2", "gene3", "gene4", "gene5", "gene6")
)

# Test loadExtInteractions
test_that("loadExtInteractions correctly outputs error", {
  expect_error(loadExtInteractions(testMOList))
  expect_error(loadExtInteractions(testMOList, NULL, NULL))
  mir <- unlist(mir2genes)
  tf <- unlist(tf2genes)
  expect_error(loadExtInteractions(testMOList, mir, tf), "list")
  mir <- list(regulator = "a")
  expect_error(loadExtInteractions(testMOList, mir))
})

test_that("loadExtInteractions load correct miRNA interactions", {
  expect_true(is.null(testMOList$extInteractions))
  testMOList <- loadExtInteractions(testMOList, mir2genes)
  expect_false(is.null(testMOList$extInteractions))
  expect_equal(testMOList$extInteractions$miR2Genes, mir2genes)
  expect_equal(testMOList$extInteractions$tf2Genes, NULL)
})

test_that("loadExtInteractions load correct TF interactions", {
  expect_true(is.null(testMOList$extInteractions))
  testMOList <- loadExtInteractions(testMOList, NULL, tf2genes)
  expect_false(is.null(testMOList$extInteractions))
  expect_equal(testMOList$extInteractions$miR2Genes, NULL)
  expect_equal(testMOList$extInteractions$tf2Genes, tf2genes)
})

test_that("loadExtInteractions load correct miRNA and TF interactions", {
  expect_true(is.null(testMOList$extInteractions))
  testMOList <- loadExtInteractions(testMOList, mir2genes, tf2genes)
  expect_false(is.null(testMOList$extInteractions))
  expect_equal(testMOList$extInteractions$miR2Genes, mir2genes)
  expect_equal(testMOList$extInteractions$tf2Genes, tf2genes)
})

test_that("loadExtInteractions coerces the correct element name", {
  expect_true(is.null(testMOList$extInteractions))
  mir2genes1 <- mir2genes
  names(mir2genes1) <- c("something", "else")
  expect_warning(testMOList <- loadExtInteractions(testMOList, mir2genes1))
  expect_false(is.null(testMOList$extInteractions))
  expect_false(any(names(testMOList$extInteractions$miR2Genes) ==
    names(mir2genes1)))
  expect_equal(
    names(testMOList$extInteractions$miR2Genes),
    c("regulator", "target")
  )
  tf2genes1 <- tf2genes
  names(tf2genes1) <- c("another", "name")
  expect_warning(testMOList <- loadExtInteractions(testMOList, NULL, tf2genes1))
  expect_false(is.null(testMOList$extInteractions))
  expect_false(any(names(testMOList$extInteractions$tf2Genes) ==
    names(tf2genes1)))
  expect_equal(
    names(testMOList$extInteractions$tf2Genes),
    c("regulator", "target")
  )
  expect_equal(names(testMOList$extInteractions$miR2Genes), NULL)
})


# Test setOmicCutoffs

test_that("setOmicCutoffs correctly outputs error", {
  expect_error(setOmicCutoffs("sdjfn"), "numeric")
  expect_silent(setOmicCutoffs())
  expect_error(setOmicCutoffs(proteomicsAdjPval = -2))
})

test_that("setOmicCutoffs correctly sets default cutoffs", {
  omiCutoffs <- setOmicCutoffs()
  expect_s3_class(omiCutoffs, "OMICutoffs")
  expect_equal(omiCutoffs$rnaAdjPval, 0.05)
  expect_equal(omiCutoffs$rnaLogFC, 0)
  expect_equal(omiCutoffs$rnaTopGenes, 1)
  expect_equal(omiCutoffs$smallRNAAdjPval, 0.05)
  expect_equal(omiCutoffs$smallRNALogFC, 0)
  expect_equal(omiCutoffs$smallRNATopGenes, 1)
  expect_equal(omiCutoffs$proteomicsAdjPval, 0.05)
  expect_equal(omiCutoffs$proteomicsLogFC, 0)
  expect_equal(omiCutoffs$atacMotifPval, NULL)
  expect_equal(omiCutoffs$atacMotifLogFC, NULL)
})

test_that("setOmicCutoffs correctly sets user-defined cutoffs", {
  omiCutoffs <- setOmicCutoffs(
    rnaAdjPval = 0.01, rnaLogFC = 1,
    rnaTopGenes = 10, smallRNAAdjPval = 0.01,
    smallRNALogFC = 1, smallRNATopGenes = 10,
    proteomicsAdjPval = 0.01, proteomicsLogFC = 1,
    atacMotifPval = 0.01,
    atacMotifLogFC = 1
  )
  expect_s3_class(omiCutoffs, "OMICutoffs")
  expect_equal(omiCutoffs$rnaAdjPval, 0.01)
  expect_equal(omiCutoffs$rnaLogFC, 1)
  expect_equal(omiCutoffs$rnaTopGenes, 10)
  expect_equal(omiCutoffs$smallRNAAdjPval, 0.01)
  expect_equal(omiCutoffs$smallRNALogFC, 1)
  expect_equal(omiCutoffs$smallRNATopGenes, 10)
  expect_equal(omiCutoffs$proteomicsAdjPval, 0.01)
  expect_equal(omiCutoffs$proteomicsLogFC, 1)
  expect_equal(omiCutoffs$atacMotifPval, 0.01)
  expect_equal(omiCutoffs$atacMotifLogFC, 1)
})

# Test constructTRN
test_that("test only RNAseq data", {
  omiCutoffs <- setOmicCutoffs()
  expect_error(constructTRN(testMOList, omiCutoffs))
})

test_that("test lacking deRNAseq data", {
  # Use example data
  data("expMOList")
  molist <- expMOList
  molist@.Data <- list()
  omiCutoffs <- setOmicCutoffs()
  expect_error(constructTRN(molist, omiCutoffs))
})

test_that("test only having deRNAseq data", {
  # Use example data
  data("expMOList")
  molist <- expMOList
  molist$DEsmallRNAseq <- NULL
  molist$annoSncRNA <- NULL
  molist$extInteractions <- NULL
  molist$DEproteomics <- NULL
  expect_error(constructTRN(molist, omiCutoffs))
  molist <- expMOList
})

test_that("test only having deRNAseq and proteomics data", {
  # Use example data
  data("expMOList")
  molist <- expMOList
  molist$DEsmallRNAseq <- NULL
  molist$annoSncRNA <- NULL
  molist$extInteractions <- NULL
  expect_error(constructTRN(molist, omiCutoffs))
})

test_that("test no external data and prediction is set false", {
  # Use example data
  data("expMOList")
  molist <- expMOList
  molist$DEproteomics <- NULL
  molist$extInteractions <- NULL
  expect_error(constructTRN(molist, omiCutoffs, predicted = FALSE))
})

test_that("test invalid target direction", {
  # Use example data
  data("expMOList")
  molist <- expMOList
  expect_error(constructTRN(molist, omiCutoffs, targetDirection = "invalid"))
})

test_that("test invalid small RNAseq type for target prediction", {
  # Use example data
  data("expMOList")
  molist <- expMOList
  expect_error(constructTRN(molist, omiCutoffs, smallRNAseqType = "invalid"))
})

test_that("test invalid small RNAseq type for target prediction", {
  # Use example data
  data("expMOList")
  molist <- expMOList
  expect_error(constructTRN(molist, omiCutoffs, smallRNAseqType = "invalid"))
})


# Perform integration test on constructTRN due to very high level of
# dependency between functions
test_that("test constructTRN", {
  # Use example data
  data("expMOList")
  omiCutoffs <- setOmicCutoffs()
  molist <- expMOList
  molist$DEsmallRNAseq <- NULL
  molist$annoSncRNA <- NULL
  molist$DEproteomics <- NULL
  molist$extInteractions <- NULL
  expect_error(constructTRN(molist, omiCutoffs))
  molist <- expMOList
  molist$DEsmallRNAseq <- NULL
  molist$annoSncRNA <- NULL
  molist$DEproteomics <- NULL
  expect_error(constructTRN(molist, omiCutoffs, predicted = FALSE))
  molist <- expMOList
  molist$DEsmallRNAseq <- NULL
  molist$annoSncRNA <- NULL
  molist$DEproteomics <- NULL
  expect_error(constructTRN(molist, omiCutoffs,
    predicted = FALSE,
    targetDirection = "both",
    smallRNAseqType = "piRNA"
  ))
  # Test with external data
  molist <- expMOList
  data("miR2Genes")
  data("tf2Genes")
  molist <- loadExtInteractions(molist, miR2Genes, tf2Genes)
  molist$DEATAC <- NULL
  expect_warning(constructTRN(molist, omiCutoffs, targetDirection = "up"))
})
