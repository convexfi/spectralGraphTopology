context("operators")
library(patrick)
library(testthat)
library(spectralGraphTopology)

# The L operator should return a symmetric positive semi-definite matrix
# with non-positive off diagonal elements and nonnegative
# diagonal elements
L_constraints <- function(Lw) {
  expect_true(isSymmetric.matrix(Lw))
  expect_true(all(diag(Lw) > 0))
  Lw_off <- Lw - diag(diag(Lw))
  expect_true(all(Lw_off <= 0))
  expect_true(all(colSums(Lw) == 0))
  expect_true(all(rowSums(Lw) == 0))
  eigen_values <- eigen(Lw, symmetric = TRUE, only.values = TRUE)
  # set to zero eigen values that are no bigger than 1e-12 in magnitude
  mask <- abs(eigen_values$values) < 1e-12
  eigen_values$values[mask] <- 0
  # verify that eigen values are nonnegative
  expect_true(all(eigen_values$values >= 0))
}


LOp <- function(w) {
  k <- length(w)
  n <- as.integer(.5 * (1 + sqrt(1 + 8 * k)))
  Lw <- matrix(0, n, n)
  for (i in ((n-1):1)) {
    j <- n - i
    Lw[i, (i+1):n] <- -tail(w[1:k], j)
    k <- k - j
  }

  Lw <- Lw + t(Lw)
  Lw <- Lw - diag(colSums(Lw))

  return(Lw)
}

LStarOp <- function(Y) {
  N <- ncol(Y)
  k <- as.integer(N * (N - 1) / 2)
  LStarY <- array(0., k)
  j <- 1
  l <- 2
  for (i in 1:k) {
    LStarY[i] <- Y[j, j] + Y[l, l] - Y[l, j] - Y[j, l]
    if (l == N) {
      j <- j + 1
      l <- j + 1
    } else {
      l <- l + 1
    }
  }
  return(LStarY)
}


LStarOpImpl <- function(Y) {
  n <- ncol(Y)
  k <- as.integer(n * (n - 1) / 2)
  LStarY <- array(0., k)
  for (i in 1:k) {
    w <- array(0., k)
    w[i] <- 1.
    Lw <- LOp(w)
    LStarY[i] <- sum(diag(t(Y) %*% Lw))
  }

  return(LStarY)
}


test_that("inverse of the L operator works", {
  w <- runif(10)
  expect_true(all(w == Linv(L(w))))
})


test_that("inverse of the A operator works", {
  w <- runif(10)
  expect_true(all(w == Ainv(A(w))))
})


with_parameters_test_that("L operator works in simple, manually verifiable cases", {
    Lw <- L(w)
    L_constraints(Lw)
    expect_true(all(Lw == answer))
  },
  cases(
        list(w = c(1, 2, 3, 4, 5, 6), answer = matrix(c(6, -1, -2, -3,
                                                       -1, 10, -4, -5,
                                                       -2, -4, 12, -6,
                                                       -3, -5, -6, 14), nrow=4)),
        list(w = c(1, 2, 3), answer = matrix(c(3, -1, -2,
                                              -1,  4, -3,
                                              -2, -3,  5), nrow=3)),
        list(w = c(1), answer = matrix(c(1, -1,
                                        -1,  1), nrow=2))
       )
)


with_parameters_test_that("A operator works in simple, manually verifiable cases", {
    expect_equal(A(w), answer)
  },
  cases(
        list(w = c(1, 2, 3, 4, 5, 6), answer = matrix(c(0, 1, 2, 3,
                                                        1, 0, 4, 5,
                                                        2, 4, 0, 6,
                                                        3, 5, 6, 0), nrow=4)),
        list(w = c(1, 2, 3), answer = matrix(c(0,  1,  2,
                                               1,  0,  3,
                                               2,  3,  0), nrow=3)),
        list(w = c(1), answer = matrix(c(0,  1,
                                         1,  0), nrow=2))
       )
)


with_parameters_test_that("L and A are linear operators", {
    w1 <- c(1, 2, 3, 4, 5, 6)
    w2 <- rev(w1)
    a <- runif(1)
    b <- runif(1)
    OPw1 <- operator(w1)
    OPw2 <- operator(w2)
    expect_equal(a * OPw1 + b * OPw2, operator(a * w1 + b * w2))
  },
  cases(list(operator = L), list(operator = A))
)


test_that("verify the Lstar operator in basic case", {
   Y <- diag(4)
   expect_true(all(LStarOp(Y) == array(2, 6)))
   expect_true(all(Lstar(Y) == array(2, 6)))
})


test_that("the composition of the Lstar and L works", {
  p <- 10
  l <- .5 * p * (p - 1)
  w <- runif(l)
  expect_true(all(Lstar(L(w)) == c(Mmat(l) %*% w)))
})


test_that("the composition of the Astar and A works", {
  p <- 10
  l <- .5 * p * (p - 1)
  w <- runif(l)
  expect_true(all(Astar(A(w)) == c(Pmat(l) %*% w)))
})


test_that("test the agreement of different implementations of the Lstar operator", {
   Y <- matrix(rnorm(16), 4, 4)
   w1 <- LStarOp(Y)
   w2 <- LStarOpImpl(Y)
   w3 <- Lstar(Y)
   expect_true(all(abs(w1 - w2) < 1e-9))
   expect_true(all(abs(w2 - w3) < 1e-9))
})


test_that("test the inner product relation between the operators L and Lstar", {
  # section 1.1 talks about an inner product equality relation
  # involving LOp and LStarOp, let's verify that
  n <- 4
  w <- c(1, 2, 3, 4, 5, 6)
  Y <- diag(n)
  Lw <- LOp(w)
  y <- LStarOp(Y)
  expect_true(sum(diag(t(Y) %*% Lw)) == w %*% y)
  Lw <- L(w)
  y <- Lstar(Y)
  expect_true(sum(diag(t(Y) %*% Lw)) == w %*% y)
})


# Neither Eigen nor base R provide the canonical vec operator,
# so we need to code our own in C++ in order to compute the matrix form of
# vec(L). Anyways, let's test our vec against matrixcalc::vec.
test_that("test our vec operator with matrixcalc::vec", {
  ncols <- sample(1:10, 1)
  nrows <- sample(1:10, 1)
  M <- matrix(runif(ncols * nrows), nrows, ncols)
  expect_true(all(vec(M) == matrixcalc::vec(M)))
})


test_that("test the composition of the operator vec and L in manually verifiable cases", {
  # n = 2
  R2 <- matrix(c(1, -1, -1, 1), 4, 1)
  expect_true(all(R2 == vecLmat(2)))
  # n = 3
  R3 <- matrix(c(1, 1, 0,
                -1, 0, 0,
                 0,-1, 0,
                -1, 0, 0,
                 1, 0, 1,
                 0, 0,-1,
                 0,-1, 0,
                 0, 0,-1,
                 0, 1, 1), 9, 3, byrow = TRUE)
  expect_true(all(R3 == vecLmat(3)))
  # n = 4
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
  expect_true(all(R4 == vecLmat(4)))
})
