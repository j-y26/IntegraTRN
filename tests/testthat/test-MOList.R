# Generate a 4x3 matrix with random integers between 1 and 10
testMatrix <- matrix(sample(1:10, 4*3, replace = TRUE), nrow = 4)
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
  expect_error(validateMOInputs(RNAseq = myMatrix, RNAGroupBy = testGroupBy,
                                smallRNAseq = NULL, smallRNAGroupBy = NULL,
                                proteomics = NULL, proteomicsGroupBy = NULL,
                                peakCond1 = NULL, peakCond2 = NULL), "NA")
  expect_error(validateMOInputs(RNAseq = testMatrix, RNAGroupBy = testGroupBy,
                                smallRNAseq = myMatrix,
                                smallRNAGroupBy = testGroupBy,
                                proteomics = NULL, proteomicsGroupBy = NULL,
                                peakCond1 = NULL, peakCond2 = NULL), "NA")
  expect_error(validateMOInputs(RNAseq = testMatrix, RNAGroupBy = testGroupBy,
                                smallRNAseq = NULL, smallRNAGroupBy = NULL,
                                proteomics = myMatrix,
                                proteomicsGroupBy = testGroupBy,
                                peakCond1 = NULL, peakCond2 = NULL), "NA")
})


