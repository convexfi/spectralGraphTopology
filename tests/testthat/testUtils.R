library(testthat)
library(spectralGraphTopology)

test_that("testBlockDiagonal", {
  N1 <- sample(1:10, 1)
  N2 <- sample(1:10, 1)
  L1 <- matrix(-1, N1, N1)
  L2 <- matrix(2, N2, N2)
  L <- rbind(cbind(L1, matrix(0, N1, N2)),
             cbind(matrix(0, N2, N1), L2))
  expect_that(all(L == blockDiag(list(L1, L2))), is_true())
})

test_that("testBlockDiagonalThrowsException", {
  N1 <- sample(1:10, 1)
  N2 <- N1 + 1
  L2 <- matrix(2, N1, N2)
  expect_error(blockDiag(list(matrix(-1, N1, N1),
                              matrix(-1, N1, N2))), "matrix 2 is not square")
})
