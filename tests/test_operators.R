library(testthat)
library(spectralGraphTopology)

test_that("test_symmetry_of_LOp", {
  # LOp should return a symmetric positive semi-definite matrix
  # with non-positive off diagonal elements and nonnegative
  # diagonal elements
  W <- LOp(w)
  expect_that(isSymmetric.matrix(W), is_true())
  expect_that(all.equal(diag(W) >= 0), is_true())
  Woff <- W - diag(diag(W))
  expect_that(all.equal(Woff <= 0), is_true())
  expect_that(all.equal(colSums(W) == 0), is_true())
  expect_that(all.equal(rowSums(W) == 0), is_true())
})
