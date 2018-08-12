library(testthat)
library(spectralGraphTopology)

LOpConstraints <- function(Lw) {
  # LOp should return a symmetric positive semi-definite matrix
  # with non-positive off diagonal elements and nonnegative
  # diagonal elements
  expect_that(isSymmetric.matrix(Lw), is_true())
  expect_that(all(diag(Lw) >= 0), is_true())
  Lw_off <- Lw - diag(diag(Lw))
  expect_that(all(Lw_off <= 0), is_true())
  expect_that(all(colSums(Lw) == 0), is_true())
  expect_that(all(rowSums(Lw) == 0), is_true())
  eigen_values <- eigen(Lw, symmetric = TRUE, only.values = TRUE)
  # set to zero eigen values that are no bigger than 1e-12 in magnitude
  mask <- abs(eigen_values$values) < 1e-12
  eigen_values$values[mask] <- 0
  # verify that eigen values are nonnegative
  expect_that(all(eigen_values$values >= 0), is_true())
}

test_that("test_LOp_order4", {
  w <- c(1, 2, 3, 4, 5, 6)
  Lw <- LOp(w, 4)
  LOpConstraints(Lw)
  expect_that(all(Lw == matrix(c(6, -1, -2, -3,
                                 -1, 10, -4, -5,
                                 -2, -4, 12, -6,
                                 -3, -5, -6, 14), nrow=4)), is_true())
})

test_that("test_LOp_order3", {
  w <- c(1, 2, 3)
  Lw <- LOp(w, 3)
  LOpConstraints(Lw)
  expect_that(all(Lw == matrix(c(3, -1, -2,
                                 -1, 4, -3,
                                 -2, -3, 5), nrow=3)), is_true())
})

test_that("test_LOp_order2", {
  w <- c(1)
  Lw <- LOp(w, 2)
  LOpConstraints(Lw)
  expect_that(all(Lw == matrix(c(1, -1, -1, 1), nrow=2)), is_true())
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

test_that("test_LStarOp", {
   # Test the LStar operator in a basic case
   Y <- diag(4)
   w <- LStarOp(Y)
   expect_that(all(w == array(2, 6)), is_true())
})

test_that("test_inner_product_relation_between_LOp_and_LStarOp", {
  # section 1.1 talks about an inner product equality relation
  # involving LOp and LStarOp, let's verify that
  n <- 4
  w <- c(1, 2, 3, 4, 5, 6)
  Lw <- LOp(w, n)
  Y <- diag(n)
  y <- LStarOp(Y)
  expect_that(sum(diag(t(Y) %*% Lw)) == w %*% y, is_true())
})
