# Purpose: Testing TRNet
# Author: Jielin Yang
# Date: 2023-11-13
# Version: 1.0
# Bugs and Issues: None

# Create a TRNet object
trnmetadata <- data.frame(
  regulator = c("A", "A", "A", "B", "B", "B"),
  target = c("C", "D", "E", "C", "D", "E"),
  regulatorType = c("TF", "TF", "TF", "miRNA", "miRNA", "miRNA")
)

# Constructor

test_that("TRNet object is created", {
  expect_equal(class(TRNet(trnmetadata, TRUE, "RNAseq"))[1], "TRNet")
})

test_that("not data.frame", {
  expect_error(TRNet(list(), TRUE, "RNAseq"))
})

test_that("not all fields in metadata", {
  df <- data.frame(
    regulator = c("A", "A", "A", "B", "B", "B"),
    target = c("C", "D", "E", "C", "D", "E")
  )
  expect_error(TRNet(df, TRUE, "RNAseq"))
})

test_that("is a list", {
  expect_true(is.list(TRNet(trnmetadata, TRUE, "RNAseq")))
  expect_true(is.list(TRNet(trnmetadata, FALSE, "RNAseq")))
  expect_true(is.list(TRNet(trnmetadata, TRUE, c("RNAseq", "ATACseq"))))
  expect_true(is.list(TRNet(trnmetadata, FALSE, c("RNAseq", "ATACseq"))))
  expect_true(is.list(TRNet(trnmetadata, TRUE, "RNAseq")@network))
})

# generatePlot

test_that("generatePlot", {
  trn <- TRNet(trnmetadata, TRUE, "RNAseq")
  trn <- generatePlot(trn)
  expect_equal(class(trn@network), "igraph")
})

test_that("correst number of edges", {
  trn <- TRNet(trnmetadata, TRUE, "RNAseq")
  trn <- generatePlot(trn)
  expect_equal(length(igraph::E(trn@network)), 6)
})

test_that("correst number of nodes", {
  trn <- TRNet(trnmetadata, TRUE, "RNAseq")
  trn <- generatePlot(trn)
  expect_equal(length(igraph::V(trn@network)), 5)
})

test_that("more omics types", {
  trn <- TRNet(trnmetadata, TRUE, c("RNAseq", "ATACseq"))
  trn <- generatePlot(trn)
  expect_equal(length(igraph::E(trn@network)), 6)
  expect_equal(length(igraph::V(trn@network)), 5)
})

# parseVertexMetadata

test_that("parse correct vertices", {
  trn <- TRNet(trnmetadata, TRUE, "RNAseq")
  trn <- generatePlot(trn)
  vertexDf <- parseVertexMetadata(trn)
  expect_equal(nrow(vertexDf), 5)
})

test_that("parse correct vertices with more omics types", {
  trn <- TRNet(trnmetadata, TRUE, c("RNAseq", "ATACseq"))
  trn <- generatePlot(trn)
  vertexDf <- parseVertexMetadata(trn)
  expect_equal(nrow(vertexDf), 5)
})

test_that("contains both regulator and target", {
  trn <- TRNet(trnmetadata, TRUE, "RNAseq")
  vertexDf <- parseVertexMetadata(trn)
  expect_true(all(c("A", "B", "C", "D", "E") %in% vertexDf$name))
})

test_that("vertex with correct types", {
  trn <- TRNet(trnmetadata, TRUE, "RNAseq")
  vertexDf <- parseVertexMetadata(trn)
  expect_equal(vertexDf[vertexDf$name == "A", "type"], "TF")
  expect_equal(vertexDf[vertexDf$name == "B", "type"], "miRNA")
  expect_equal(vertexDf[vertexDf$name == "C", "type"], "gene")
  expect_equal(vertexDf[vertexDf$name == "D", "type"], "gene")
  expect_equal(vertexDf[vertexDf$name == "E", "type"], "gene")
})

# plotNetwork

test_that("plotNetwork static", {
  trn <- TRNet(trnmetadata, TRUE, "RNAseq")
  network <- plotNetwork(trn)
  expect_null(network)
})

test_that("plotNetwork interactive", {
  trn <- TRNet(trnmetadata, TRUE, "RNAseq")
  network <- plotNetwork(trn, interactive = TRUE)
  expect_equal(class(network), c("forceNetwork", "htmlwidget"))
})

# exportIgraph

test_that("exportIgraph", {
  trn <- TRNet(trnmetadata, TRUE, "RNAseq")
  igraph <- exportIgraph(trn)
  expect_equal(class(igraph), "igraph")
})

# exportEdgeSet

test_that("exportEdgeSet", {
  trn <- TRNet(trnmetadata, TRUE, "RNAseq")
  edgeSet <- exportEdgeSet(trn)
  expect_equal(class(edgeSet), "data.frame")
  expect_equal(nrow(edgeSet), 6)
  expect_equal(ncol(edgeSet), 3)
})

# [END]
