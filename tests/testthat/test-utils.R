context("utility functions")
library(testthat)
library(spectralGraphTopology)


test_that("consistency of blockDiag function", {
  N1 <- sample(1:10, 1)
  N2 <- sample(1:10, 1)
  L1 <- matrix(-1, N1, N1)
  L2 <- matrix(2, N2, N2)
  L <- rbind(cbind(L1, matrix(0, N1, N2)),
             cbind(matrix(0, N2, N1), L2))
  expect_true(all(L == block_diag(L1, L2)))
})


test_that("upper_view_vec works", {
  a <- runif(1)
  M <- a * matrix(c(6, -1, -2, -3,
                    -1, 10, -4, -5,
                    -2, -4, 12, -6,
                    -3, -5, -6, 14))
  v <- upper_view_vec(M)
  expect_true(all(v == - a * c(1:6)))
})
