library(testthat)
library(spectralGraphTopology)

LOpConstraints <- function(Lw) {
  # The L operator should return a symmetric positive semi-definite matrix
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

test_that("test_Linv", {
  w <- runif(10)
  Lw <- L(w)
  w_test <- Linv(Lw)
  expect_that(all(w == w_test), is_true())
})

test_that("test_LOp_order4", {
  w <- c(1, 2, 3, 4, 5, 6)
  answer <- matrix(c(6, -1, -2, -3,
                     -1, 10, -4, -5,
                     -2, -4, 12, -6,
                     -3, -5, -6, 14), nrow=4)
  Lw <- LOp(w)
  LOpConstraints(Lw)
  expect_that(all(Lw == answer), is_true())

  Lw <- L(w)
  LOpConstraints(Lw)
  expect_that(all(Lw == answer), is_true())
})

test_that("test_LOp_order3", {
  w <- c(1, 2, 3)
  answer <- matrix(c(3, -1, -2,
                     -1, 4, -3,
                     -2, -3, 5), nrow=3)
  Lw <- LOp(w)
  LOpConstraints(Lw)
  expect_that(all(Lw == answer), is_true())

  Lw <- L(w)
  LOpConstraints(Lw)
  expect_that(all(Lw == answer), is_true())
})

test_that("test_LOp_order2", {
  w <- c(1)
  answer <- matrix(c(1, -1, -1, 1), nrow=2)
  Lw <- LOp(w)
  LOpConstraints(Lw)
  expect_that(all(Lw == answer), is_true())
  Lw <- L(w)
  LOpConstraints(Lw)
  expect_that(all(Lw == answer), is_true())
})

test_that("test_linearity_of_LOp", {
  # Verify that the implemented L operator is indeed linear
  w1 <- c(1, 2, 3, 4, 5, 6)
  w2 <- rev(w1)
  a <- runif(1)
  b <- runif(1)
  Lw1 <- LOp(w1)
  Lw2 <- LOp(w2)
  expect_that(all((a * Lw1 + b * Lw2) == LOp(a * w1 + b * w2)), is_true())
  Lw1 <- L(w1)
  Lw2 <- L(w2)
  expect_that(all((a * Lw1 + b * Lw2) == L(a * w1 + b * w2)), is_true())
})

test_that("test_LStarOp", {
   # Test the LStar operator in a basic case
   Y <- diag(4)
   w <- LStarOp(Y)
   expect_that(all(w == array(2, 6)), is_true())

   w <- Lstar(Y)
   expect_that(all(w == array(2, 6)), is_true())
})

test_that("test_LStarOp_random", {
   # Test the LStar operator implementations are compatible among themselves
   Y <- matrix(rnorm(16), 4, 4)
   w1 <- LStarOp(Y)
   w2 <- LStarOpImpl(Y)
   w3 <- Lstar(Y)
   expect_that(all(abs(w1 - w2) < 1e-6), is_true())
   expect_that(all(abs(w2 - w3) < 1e-6), is_true())
})

test_that("test_inner_product_relation_between_LOp_and_LStarOp", {
  # section 1.1 talks about an inner product equality relation
  # involving LOp and LStarOp, let's verify that
  n <- 4
  w <- c(1, 2, 3, 4, 5, 6)
  Y <- diag(n)
  Lw <- LOp(w)
  y <- LStarOp(Y)
  expect_that(sum(diag(t(Y) %*% Lw)) == w %*% y, is_true())

  Lw <- L(w)
  y <- Lstar(Y)
  expect_that(sum(diag(t(Y) %*% Lw)) == w %*% y, is_true())
})

# Unfortunately, neither Eigen nor base R provide the canonical vec operator,
# so we need to code our own in C++ in order to compute the matrix form of
# vec(L). Anyways, let's test our vec against matrixcalc::vec.
test_that("test_vec_operator_against_matrixcalc", {
    ncols <- sample(1:10, 1)
    nrows <- sample(1:10, 1)
    M <- matrix(runif(ncols * nrows), nrows, ncols)
    expect_that(all(vec(M) == matrixcalc::vec(M)), is_true())
})

test_that("test_vec(L)_operator", {
    # N = 2
    R2 <- matrix(c(1, -1, -1, 1), 4, 1)
    expect_that(all(R2 == vecLmat(2)), is_true())

    # N = 3
    R3 <- matrix(c(1, 1, 0,
                  -1, 0, 0,
                   0,-1, 0,
                  -1, 0, 0,
                   1, 0, 1,
                   0, 0,-1,
                   0,-1, 0,
                   0, 0,-1,
                   0, 1, 1), 9, 3, byrow = TRUE)
    expect_that(all(R3 == vecLmat(3)), is_true())

    # N = 4
    R4 <- matrix(c(1, 1, 1, 0, 0, 0,
                  -1, 0, 0, 0, 0, 0,
                   0,-1, 0, 0, 0, 0,
                   0, 0,-1, 0, 0, 0,
                  -1, 0, 0, 0, 0, 0,
                   1, 0, 0, 1, 1, 0,
                   0, 0, 0,-1, 0, 0,
                   0, 0, 0, 0,-1, 0,
                   0,-1, 0, 0, 0, 0,
                   0, 0, 0,-1, 0, 0,
                   0, 1, 0, 1, 0, 1,
                   0, 0, 0, 0, 0,-1,
                   0, 0,-1, 0, 0, 0,
                   0, 0, 0, 0,-1, 0,
                   0, 0, 0, 0, 0,-1,
                   0, 0, 1, 0, 1, 1), 16, 6, byrow = TRUE)
    expect_that(all(R4 == vecLmat(4)), is_true())
})
