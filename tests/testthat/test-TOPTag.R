# Purpose: Testing TOPTag
# Author: Jielin Yang
# Date: 2023-11-13
# Version: 1.0
# Bugs and Issues: None


# Define usable objects
deResult <- data.frame(
  log2FoldChange = c(1, 2, 3, -4, -5),
  pvalue = c(0.0001, 0.0002, 0.0003, 0.0004, 0.0005),
  padj = c(0.0004, 0.0005, 0.0006, 0.0007, 0.0008)
)

deResult2 <- data.frame(
  logFC = c(1, 2, 3, -4, -5),
  PValue = c(0.0001, 0.0002, 0.0003, 0.0004, 0.0005),
  FDR = c(0.0004, 0.0005, 0.0006, 0.0007, 0.0008)
)

genes <- c("gene1", "gene2", "gene3", "gene4", "gene5")

rownames(deResult) <- genes
rownames(deResult2) <- genes
normCounts <- data.frame(
  sample1 = c(1, 2, 3, 4, 5),
  sample2 = c(6, 7, 8, 9, 10),
  sample3 = c(11, 12, 13, 14, 15),
  sample4 = c(16, 17, 18, 19, 20),
  sample5 = c(21, 22, 23, 24, 25)
)
normCounts <- normCounts / 1000
normCounts <- as.matrix(normCounts)
rownames(normCounts) <- genes

deTag <- DETag(deResult, "DESeq2", normCounts)
deTag2 <- DETag(deResult2, "edgeR", normCounts)

# Constructor

test_that("do not filter deseq2", {
  toptag <- TOPTag(deTag,
    logFCCutoff = 0, pCutoff = 1,
    topGenes = 1, direction = "both"
  )
  expect_equal(toptag@topGenes, 1)
  expect_equal(toptag@logFCCutoff, 0)
  expect_equal(toptag@pCutoff, 1)
  deRes <- deResult
  colnames(deRes) <- COUNT_DEFIELDS
  for (gene in genes) {
    expect_equal(toptag@DEResult[gene, 1:3], deRes[gene, ])
  }
  expect_true(all(COUNT_DEFIELDS %in% colnames(toptag@DEResult)))
  expect_true(all(c("piValue", "rank") %in% colnames(toptag@DEResult)))
  for (gene in genes) {
    expect_equal(
      -log10(toptag@DEResult$padj[gene]) *
        sign(toptag@DEResult$logFC[gene]),
      toptag@DEResult$piValue[gene]
    )
  }
  expect_equal(sum(toptag@DEResult$rank < 0), 2)
  expect_equal(sum(toptag@DEResult$rank > 0), 3)
  for (gene in genes) {
    expect_equal(toptag@normalizedCounts[gene, ], normCounts[gene, ])
  }
  expect_equal(toptag@method, "DESeq2")
})

test_that("do not filter edgeR", {
  toptag <- TOPTag(deTag2,
    logFCCutoff = 0, pCutoff = 1,
    topGenes = 1, direction = "both"
  )
  expect_equal(toptag@topGenes, 1)
  expect_equal(toptag@logFCCutoff, 0)
  expect_equal(toptag@pCutoff, 1)
  deRes <- deResult2
  colnames(deRes) <- COUNT_DEFIELDS
  for (gene in genes) {
    expect_equal(toptag@DEResult[gene, 1:3], deRes[gene, ])
  }
  expect_true(all(COUNT_DEFIELDS %in% colnames(toptag@DEResult)))
  expect_true(all(c("piValue", "rank") %in% colnames(toptag@DEResult)))
  for (gene in genes) {
    expect_equal(
      -log10(toptag@DEResult$padj[gene]) *
        sign(toptag@DEResult$logFC[gene]),
      toptag@DEResult$piValue[gene]
    )
  }
  expect_equal(sum(toptag@DEResult$rank < 0), 2)
  expect_equal(sum(toptag@DEResult$rank > 0), 3)
  for (gene in genes) {
    expect_equal(toptag@normalizedCounts[gene, ], normCounts[gene, ])
  }
  expect_equal(toptag@method, "edgeR")
})

test_that("filter deseq2", {
  toptag <- TOPTag(deTag,
    logFCCutoff = 1, pCutoff = 0.0007,
    topGenes = 1, direction = "both"
  )
  expect_equal(toptag@topGenes, 1)
  expect_equal(toptag@logFCCutoff, 1)
  expect_equal(toptag@pCutoff, 0.0007)
  deRes <- deResult
  colnames(deRes) <- COUNT_DEFIELDS
  expect_equal(nrow(toptag@DEResult), 2)
  expect_true(all(rownames(toptag@DEResult) %in% c("gene2", "gene3")))
  expect_true(all(COUNT_DEFIELDS %in% colnames(toptag@DEResult)))
  expect_true(all(c("piValue", "rank") %in% colnames(toptag@DEResult)))
  expect_equal(sum(toptag@DEResult$rank < 0), 0)
  expect_equal(sum(toptag@DEResult$rank > 0), 2)
  expect_equal(toptag@method, "DESeq2")
})

test_that("filter edgeR", {
  toptag <- TOPTag(deTag2,
    logFCCutoff = 1, pCutoff = 0.0007,
    topGenes = 1, direction = "both"
  )
  expect_equal(toptag@topGenes, 1)
  expect_equal(toptag@logFCCutoff, 1)
  expect_equal(toptag@pCutoff, 0.0007)
  deRes <- deResult2
  colnames(deRes) <- COUNT_DEFIELDS
  expect_equal(nrow(toptag@DEResult), 2)
  expect_true(all(rownames(toptag@DEResult) %in% c("gene2", "gene3")))
  expect_true(all(COUNT_DEFIELDS %in% colnames(toptag@DEResult)))
  expect_true(all(c("piValue", "rank") %in% colnames(toptag@DEResult)))
  expect_equal(sum(toptag@DEResult$rank < 0), 0)
  expect_equal(sum(toptag@DEResult$rank > 0), 2)
  expect_equal(toptag@method, "edgeR")
})

test_that("wrong method", {
  deTag3 <- DETag(deResult, "DESeq2", normCounts)
  deTag3@method <- "edgeR" # force it to be a wrong method
  expect_error(TOPTag(deTag3,
    logFCCutoff = 1, pCutoff = 0.0007,
    topGenes = 1, direction = "both"
  ))
})

test_that("wrong direction", {
  expect_error(TOPTag(deTag,
    logFCCutoff = 1, pCutoff = 0.0007,
    topGenes = 1, direction = "wrong"
  ))
})

test_that("test topGenes is a fraction", {
  toptag <- TOPTag(deTag,
    logFCCutoff = 0, pCutoff = 1,
    topGenes = 0.5, direction = "both"
  )
  expect_equal(toptag@topGenes, 0.5)
  expect_equal(toptag@logFCCutoff, 0)
  expect_equal(toptag@pCutoff, 1)
  deRes <- deResult
  colnames(deRes) <- COUNT_DEFIELDS
  expect_equal(nrow(toptag@DEResult), 3)
  expect_true(all(rownames(toptag@DEResult) %in% c("gene1", "gene2", "gene5")))
  expect_true(all(COUNT_DEFIELDS %in% colnames(toptag@DEResult)))
  expect_true(all(c("piValue", "rank") %in% colnames(toptag@DEResult)))
  expect_equal(sum(toptag@DEResult$rank < 0), 1)
  expect_equal(sum(toptag@DEResult$rank > 0), 2)
  expect_equal(toptag@method, "DESeq2")
})

test_that("test topGenes is a number", {
  toptag <- TOPTag(deTag,
    logFCCutoff = 0, pCutoff = 1,
    topGenes = 2, direction = "both"
  )
  expect_equal(toptag@topGenes, 2)
  expect_equal(toptag@logFCCutoff, 0)
  expect_equal(toptag@pCutoff, 1)
  deRes <- deResult
  colnames(deRes) <- COUNT_DEFIELDS
  expect_equal(nrow(toptag@DEResult), 2)
  expect_true(all(rownames(toptag@DEResult) %in% c("gene1", "gene5")))
  expect_true(all(COUNT_DEFIELDS %in% colnames(toptag@DEResult)))
  expect_true(all(c("piValue", "rank") %in% colnames(toptag@DEResult)))
  expect_equal(sum(toptag@DEResult$rank < 0), 1)
  expect_equal(sum(toptag@DEResult$rank > 0), 1)
  expect_equal(toptag@method, "DESeq2")
})


# filterGenes

test_that("filter genes, all genes", {
  toptag <- TOPTag(deTag,
    logFCCutoff = 0, pCutoff = 1,
    topGenes = 1, direction = "both"
  )
  expect_equal(toptag@topGenes, 1)
  expect_equal(toptag@logFCCutoff, 0)
  expect_equal(toptag@pCutoff, 1)
  expect_equal(nrow(toptag@DEResult), 5)
  filteredToptag <- filterGenes(toptag, genes)
  expect_equal(nrow(filteredToptag@DEResult), 5)
  expect_true(all(rownames(filteredToptag@DEResult) %in% genes))
  expect_equal(filteredToptag@topGenes, 1)
  expect_equal(filteredToptag@logFCCutoff, 0)
  expect_equal(filteredToptag@pCutoff, 1)
  expect_equal(filteredToptag@method, "DESeq2")
})

test_that("filter genes, some genes", {
  toptag <- TOPTag(deTag,
    logFCCutoff = 0, pCutoff = 1,
    topGenes = 1, direction = "both"
  )
  expect_equal(toptag@topGenes, 1)
  expect_equal(toptag@logFCCutoff, 0)
  expect_equal(toptag@pCutoff, 1)
  expect_equal(nrow(toptag@DEResult), 5)
  filteredToptag <- filterGenes(toptag, c("gene1", "gene2", "gene3"))
  expect_equal(nrow(filteredToptag@DEResult), 3)
  expect_true(all(rownames(filteredToptag@DEResult) %in%
    c("gene1", "gene2", "gene3")))
  expect_equal(filteredToptag@topGenes, 1)
  expect_equal(filteredToptag@logFCCutoff, 0)
  expect_equal(filteredToptag@pCutoff, 1)
  expect_equal(filteredToptag@method, "DESeq2")
})

test_that("filter genes, no genes", {
  toptag <- TOPTag(deTag,
    logFCCutoff = 0, pCutoff = 1,
    topGenes = 1, direction = "both"
  )
  expect_equal(toptag@topGenes, 1)
  expect_equal(toptag@logFCCutoff, 0)
  expect_equal(toptag@pCutoff, 1)
  expect_equal(nrow(toptag@DEResult), 5)
  filteredToptag <- filterGenes(toptag, c())
  expect_equal(nrow(filteredToptag@DEResult), 0)
  expect_equal(filteredToptag@topGenes, 1)
  expect_equal(filteredToptag@logFCCutoff, 0)
  expect_equal(filteredToptag@pCutoff, 1)
  expect_equal(filteredToptag@method, "DESeq2")
})

# exportDE

test_that("get correct columns", {
  toptag <- TOPTag(deTag,
    logFCCutoff = 0, pCutoff = 1,
    topGenes = 1, direction = "both"
  )
  expect_equal(nrow(exportDE(toptag)), 5)
  expect_true(all(rownames(exportDE(toptag)) %in% genes))
  expect_true(all(c("logFC", "pvalue", "padj", "piValue", "rank") %in%
    colnames(exportDE(toptag))))
})

# [END]
