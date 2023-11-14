# Purpose: Testing PEAKTag
# Author: Jielin Yang
# Date: 2023-11-13
# Version: 1.0
# Bugs and Issues: None

# Generate some simulated peak data
deResult <- data.frame(
  chr = c("chr1", "chr2", "chr3", "chr4", "chr5"),
  start = c(1, 2, 3, 4, 5),
  end = c(2, 3, 4, 5, 6),
  Condition = c("-", "-", "+", "+", "+")
)

# Create a DETag for peaks
detag <- DETag(deResult, "GRangesATAC")

# Constructor

test_that("constructor", {
  expect_s4_class(PEAKTag(detag), "PEAKTag")
})

test_that("construct with wrong type of data", {
  detag2 <- detag
  detag2@method <- "something wrong"
  expect_error(PEAKTag(detag2))
})

test_that("wrong TxDB", {
  expect_error(PEAKTag(detag, TxDB = "something wrong"))
})

test_that("wrong annoDB", {
  expect_error(PEAKTag(detag, annoDB = org.Hs.eg.db::org.Hs.eg.db))
})


# convert to GRanges

test_that("convert to GRanges", {
  peaktag <- PEAKTag(detag)
  expect_s4_class(asGRanges(peaktag), "GRanges")
})


# convert to data frame

test_that("convert to data frame", {
  peaktag <- PEAKTag(detag)
  expect_true(is.data.frame(as.data.frame(peaktag)))
})
