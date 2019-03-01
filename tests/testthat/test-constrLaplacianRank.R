context("Jordan's constrained Laplacian rank")
library(testthat)
library(spectralGraphTopology)

test_that("pairwise row norm", {
  n <- sample(c(3:10), size = 1)
  M <- matrix(runif(n ^ 2), n, n)
  V <- pairwise_matrix_rownorm(M)
  expect_that(isSymmetric(V), is_true())
  expect_equal(diag(V), rep(0, n))
})

test_that("Jordan's initial graph", {
  m <- sample(c(10:20), size = 1)
  n <- sample(c(10:20), size = 1)
  Y <- matrix(runif(n * m), n, m)
  A <- initial_graph(Y, m = 5)
  expect_that(isSymmetric(A), is_true())
  expect_equal(diag(A), rep(0, n))
  expect_equal(rowSums(A), rep(1, n))
})

