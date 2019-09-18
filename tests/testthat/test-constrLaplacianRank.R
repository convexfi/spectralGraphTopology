context("Jordan's constrained Laplacian rank")
library(testthat)
library(spectralGraphTopology)

test_that("pairwise row norm", {
  n <- sample(c(3:10), size = 1)
  m <- sample(c(3:10), size = 1)
  M <- matrix(runif(n * m), n, m)
  V <- pairwise_matrix_rownorm(M)
  expect_equal(c(n, n), dim(V))
  expect_true(isSymmetric(V))
  expect_equal(diag(V), rep(0, n))

  M <- rbind(c(1, 2, 3),
             c(3, 2, 0),
             c(1, 1, 1))
  V <- pairwise_matrix_rownorm(M)
  ANS <- rbind(c(0, 13, 5),
              c(13,  0, 6),
               c(5,  6, 0))
  expect_equal(V, ANS)
})

test_that("Jordan's initial graph", {
  m <- sample(c(10:20), size = 1)
  n <- 5
  Y <- matrix(runif(n * m), n, m)
  A <- build_initial_graph(Y, m = 2)
  expect_equal(c(n, n), dim(A))
  expect_equal(diag(A), rep(0, n))
  expect_equal(rowSums(A), rep(1, n))
})

