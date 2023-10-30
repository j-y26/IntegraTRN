# Generate a 4x3 matrix with random integers between 1 and 10
testMatrix <- matrix(sample(1:10, 4 * 3, replace = TRUE), nrow = 4)
colnames(testMatrix) <- c("sample1", "sample2", "sample3")
rownames(testMatrix) <- paste0("gene_", seq_len(4))
testGroupBy <- c("A", "B", "C")

# Testing validateMatrix helper function
test_that("validateMatrix throws error if matrix contains NA", {
  myMatrix <- testMatrix
  myMatrix[1, 1] <- NA
  expect_error(validateMatrix(myMatrix, testGroupBy), "NA")
})

test_that("validateMatrix throws error if matrix contains non-numeric values", {
  myMatrix <- testMatrix
  myMatrix[1, 1] <- "A"
  expect_error(validateMatrix(myMatrix, testGroupBy), "numeric")
})

test_that("validateMatrix throws error if number of columns is not equal to
  length of groupBy", {
  myMatrix <- testMatrix
  myMatrix <- myMatrix[, 1:2]
  expect_error(validateMatrix(myMatrix, testGroupBy), "match")
})

test_that("validateMatrix throws error if groupBy contains NA", {
  myGroupBy <- testGroupBy
  myGroupBy[1] <- NA
  expect_error(validateMatrix(testMatrix, myGroupBy), "provide")
  myGroupBy <- NA
  expect_error(validateMatrix(testMatrix, myGroupBy), "provide")
})

# Testing validateMOInputs helper function
test_that("validateMOInputs throws error if matrix contains NA", {
  myMatrix <- testMatrix
  myMatrix[1, 1] <- NA
  expect_error(validateMOInputs(
    RNAseq = myMatrix, RNAGroupBy = testGroupBy,
    smallRNAseq = NULL, smallRNAGroupBy = NULL,
    proteomics = NULL, proteomicsGroupBy = NULL,
    peakCond1 = NULL, peakCond2 = NULL
  ), "NA")
  expect_error(validateMOInputs(
    RNAseq = testMatrix, RNAGroupBy = testGroupBy,
    smallRNAseq = myMatrix,
    smallRNAGroupBy = testGroupBy,
    proteomics = NULL, proteomicsGroupBy = NULL,
    peakCond1 = NULL, peakCond2 = NULL
  ), "NA")
  expect_error(validateMOInputs(
    RNAseq = testMatrix, RNAGroupBy = testGroupBy,
    smallRNAseq = NULL, smallRNAGroupBy = NULL,
    proteomics = myMatrix,
    proteomicsGroupBy = testGroupBy,
    peakCond1 = NULL, peakCond2 = NULL
  ), "NA")
})

test_that("validateMOInputs throws error if matrix contains non-numeric values", {
  myMatrix <- testMatrix
  myMatrix[1, 1] <- "A"
  expect_error(
    validateMOInputs(
      RNAseq = myMatrix,
      RNAGroupBy = testGroupBy,
      smallRNAseq = NULL,
      smallRNAGroupBy = NULL,
      proteomics = NULL,
      proteomicsGroupBy = NULL,
      peakCond1 = NULL, peakCond2 = NULL
    ),
    "numeric"
  )
  expect_error(
    validateMOInputs(
      RNAseq = testMatrix,
      RNAGroupBy = testGroupBy,
      smallRNAseq = myMatrix,
      smallRNAGroupBy = testGroupBy,
      proteomics = NULL,
      proteomicsGroupBy = NULL,
      peakCond1 = NULL, peakCond2 = NULL
    ),
    "numeric"
  )
  expect_error(
    validateMOInputs(
      RNAseq = testMatrix,
      RNAGroupBy = testGroupBy,
      smallRNAseq = NULL,
      smallRNAGroupBy = NULL,
      proteomics = myMatrix,
      proteomicsGroupBy = testGroupBy,
      peakCond1 = NULL, peakCond2 = NULL
    ),
    "numeric"
  )
})

test_that("validateMOInputs throws error if number of columns is not equal to
  length of groupBy", {
  myMatrix <- testMatrix
  myMatrix <- myMatrix[, 1:2]
  expect_error(validateMOInputs(
    RNAseq = myMatrix, RNAGroupBy = testGroupBy,
    smallRNAseq = NULL, smallRNAGroupBy = NULL,
    proteomics = NULL, proteomicsGroupBy = NULL,
    peakCond1 = NULL, peakCond2 = NULL
  ), "match")
  expect_error(validateMOInputs(
    RNAseq = testMatrix, RNAGroupBy = testGroupBy,
    smallRNAseq = myMatrix,
    smallRNAGroupBy = testGroupBy,
    proteomics = NULL, proteomicsGroupBy = NULL,
    peakCond1 = NULL, peakCond2 = NULL
  ), "match")
  expect_error(validateMOInputs(
    RNAseq = testMatrix, RNAGroupBy = testGroupBy,
    smallRNAseq = NULL, smallRNAGroupBy = NULL,
    proteomics = myMatrix,
    proteomicsGroupBy = testGroupBy,
    peakCond1 = NULL, peakCond2 = NULL
  ), "match")
})

test_that("validateMOInputs throws error if groupBy contains NA", {
  myGroupBy <- testGroupBy
  myGroupBy[1] <- NA
  expect_error(validateMOInputs(
    RNAseq = testMatrix, RNAGroupBy = myGroupBy,
    smallRNAseq = NULL, smallRNAGroupBy = NULL,
    proteomics = NULL, proteomicsGroupBy = NULL,
    peakCond1 = NULL, peakCond2 = NULL
  ), "provide")
  expect_error(validateMOInputs(
    RNAseq = testMatrix, RNAGroupBy = testGroupBy,
    smallRNAseq = testMatrix,
    smallRNAGroupBy = myGroupBy,
    proteomics = NULL, proteomicsGroupBy = NULL,
    peakCond1 = NULL, peakCond2 = NULL
  ), "provide")
  expect_error(validateMOInputs(
    RNAseq = testMatrix, RNAGroupBy = testGroupBy,
    smallRNAseq = NULL, smallRNAGroupBy = NULL,
    proteomics = testMatrix,
    proteomicsGroupBy = myGroupBy,
    peakCond1 = NULL, peakCond2 = NULL
  ), "provide")
  expect_error(validateMOInputs(
    RNAseq = testMatrix, RNAGroupBy = testGroupBy,
    smallRNAseq = NULL, smallRNAGroupBy = NULL,
    proteomics = testMatrix,
    proteomicsGroupBy = NULL,
    peakCond1 = NULL, peakCond2 = NULL
  ), "match")
})

test_that("validateMOInputs throws error if some data does not contain
  replicate values", {
  myMatrix <- testMatrix
  myMatrix <- matrix(testMatrix[, 1], ncol = 1)
  expect_error(
    validateMOInputs(
      RNAseq = myMatrix, RNAGroupBy = 1,
      smallRNAseq = NULL, smallRNAGroupBy = NULL,
      proteomics = NULL, proteomicsGroupBy = NULL,
      peakCond1 = NULL, peakCond2 = NULL
    ),
    "replicates"
  )
  expect_error(
    validateMOInputs(
      RNAseq = testMatrix, RNAGroupBy = testGroupBy,
      smallRNAseq = myMatrix,
      smallRNAGroupBy = 1,
      proteomics = NULL, proteomicsGroupBy = NULL,
      peakCond1 = NULL, peakCond2 = NULL
    ),
    "replicates"
  )
  expect_error(
    validateMOInputs(
      RNAseq = testMatrix, RNAGroupBy = testGroupBy,
      smallRNAseq = NULL, smallRNAGroupBy = NULL,
      proteomics = myMatrix,
      proteomicsGroupBy = 1,
      peakCond1 = NULL, peakCond2 = NULL
    ),
    "replicates"
  )
})

# Testing MOList object creation
test_that("MOList object initial creation with only RNAseq data", {
  myMOList <- MOList(RNAseq = testMatrix, RNAGroupBy = testGroupBy)
  expect_equal(myMOList@RNAseq, testMatrix)
  expect_equal(myMOList@RNAseqSamples$groupBy, testGroupBy)
  expect_equal(length(myMOList@RNAseqSamples$samples), ncol(myMOList@RNAseq))
  expect_true(is.matrix(myMOList@smallRNAseq))
  expect_equal(myMOList@smallRNAseqSamples$samples, NULL)
  expect_equal(myMOList@smallRNAseqSamples$groupBy, NULL)
  expect_true(is.matrix(myMOList@proteomics))
  expect_equal(myMOList@proteomicsSamples$samples, NULL)
  expect_equal(myMOList@proteomicsSamples$groupBy, NULL)
  expect_equal(myMOList@ATACpeaks$peaksCond1, NULL)
  expect_equal(myMOList@ATACpeaks$peaksCond2, NULL)
})

test_that("MOList object creation without RNAseq data", {
  expect_error(MOList(), "either")
})

test_that("MOList object creation with additional data", {
  myMOList <- MOList(
    RNAseq = testMatrix, RNAGroupBy = testGroupBy,
    smallRNAseq = testMatrix, smallRNAGroupBy = testGroupBy,
    proteomics = testMatrix, proteomicsGroupBy = testGroupBy
  )
  expect_equal(myMOList@RNAseq, testMatrix)
  expect_equal(myMOList@RNAseqSamples$groupBy, testGroupBy)
  expect_equal(length(myMOList@RNAseqSamples$samples), ncol(myMOList@RNAseq))
  expect_equal(myMOList@smallRNAseq, testMatrix)
  expect_equal(myMOList@smallRNAseqSamples$groupBy, testGroupBy)
  expect_equal(
    length(myMOList@smallRNAseqSamples$samples),
    ncol(myMOList@smallRNAseq)
  )
  expect_equal(myMOList@proteomics, testMatrix)
  expect_equal(myMOList@proteomicsSamples$groupBy, testGroupBy)
  expect_equal(
    length(myMOList@proteomicsSamples$samples),
    ncol(myMOList@proteomics)
  )
})

test_that("Appending data to MOList object", {
  myMOList <- MOList(RNAseq = testMatrix, RNAGroupBy = testGroupBy)
  expect_equal(myMOList@RNAseq, testMatrix)
  expect_equal(myMOList@RNAseqSamples$groupBy, testGroupBy)
  expect_equal(length(myMOList@RNAseqSamples$samples), ncol(myMOList@RNAseq))
  expect_true(ncol(myMOList@smallRNAseq) == 0)
  myMOList <- MOList(myMOList,
    smallRNAseq = testMatrix,
    smallRNAGroupBy = testGroupBy
  )
  expect_equal(myMOList@smallRNAseq, testMatrix)
  expect_equal(myMOList@smallRNAseqSamples$groupBy, testGroupBy)
  expect_equal(
    length(myMOList@smallRNAseqSamples$samples),
    ncol(myMOList@smallRNAseq)
  )
  expect_true(ncol(myMOList@proteomics) == 0)
  myMOList <- MOList(myMOList,
    proteomics = testMatrix,
    proteomicsGroupBy = testGroupBy
  )
  expect_equal(myMOList@proteomics, testMatrix)
  expect_equal(myMOList@proteomicsSamples$groupBy, testGroupBy)
  expect_equal(
    length(myMOList@proteomicsSamples$samples),
    ncol(myMOList@proteomics)
  )
  expect_false(is.data.frame(myMOList@ATACpeaks$peaksCond1))
  expect_false(is.data.frame(myMOList@ATACpeaks$peaksCond2))
})

# Testing MOList object validation functions
test_that("MOList object validation functions", {
  myMOList <- MOList(RNAseq = testMatrix, RNAGroupBy = testGroupBy,
                     smallRNAseq = testMatrix, smallRNAGroupBy = testGroupBy,
                     proteomics = testMatrix, proteomicsGroupBy = testGroupBy
  )
  # Temper the object illegally
  testMOList <- myMOList
  testMOList@RNAseq[1, 1] <- NA
  expect_error(validateMOList(testMOList), "NA")
  testMOList <- myMOList
  testMOList@RNAseq[1, 1] <- "A"
  expect_error(validateMOList(testMOList), "numeric")
  testMOList <- myMOList
  testMOList@RNAseq <- testMatrix[, 1:2]
  expect_error(validateMOList(testMOList), "match")
  testMOList <- myMOList
  testMOList@RNAseqSamples$groupBy[1] <- NA
  expect_error(validateMOList(testMOList), "provide")
  testMOList <- myMOList
  testMOList@RNAseqSamples$groupBy <- testMOList@RNAseqSamples$groupBy[-1]
  expect_error(validateMOList(testMOList), "correct")
})
