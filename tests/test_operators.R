library(testthat)
library(spectralGraphTopology)

test_that("test_symmetry_of_LOp", {
  # LOp should return a symmetric positive semi-definite matrix
  # with non-positive off diagonal elements and nonnegative
  # diagonal elements
  w <- c(1, 2, 3, 4, 5, 6)
  Lw <- LOp(w, 4)
  expect_that(isSymmetric.matrix(Lw), is_true())
  expect_that(all(diag(Lw) >= 0), is_true())
  Lw_off <- Lw - diag(diag(Lw))
  expect_that(all(Lw_off <= 0), is_true())
  expect_that(all(colSums(Lw) == 0), is_true())
  expect_that(all(rowSums(Lw) == 0), is_true())
})

test_that("test_linearity_of_LOp", {
  # Verify that the implemented L operator is indeed linear
  w1 <- c(1, 2, 3, 4, 5, 6)
  w2 <- rev(w1)
  a <- runif(1)
  b <- runif(1)
  Lw1 <- LOp(w1, 4)
  Lw2 <- LOp(w2, 4)
  expect_that(all((a * Lw1 + b * Lw2) == LOp(a * w1 + b * w2, 4)), is_true())
})

test_that("test_inner_product_relation_between_LOp_and_LStarOp", {
  # section 1.1 talks about a inner product equality relation
  # involving LOp and LStarOp, let's verify that
  Lw <- LOp(w)
  y <- LStarOp(Y)
  expect_that(sum(diag(t(Y) %*% Lw)) == w %*% y, is_true())
})
