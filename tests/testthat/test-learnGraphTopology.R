context("Laplacian Matrix Estimation")
library(testthat)
library(spectralGraphTopology)


test_that("learn_laplacian_matrix with single component random graph", {
  w <- sample(1:10, 6)
  Lw <- L(w)
  n <- ncol(Lw)
  Y <- MASS::mvrnorm(1000, rep(0, n), MASS::ginv(Lw))
  res <- learn_laplacian_matrix(cov(Y))
  expect_that(relativeError(Lw, res$Lw) < 1e-1, is_true())
  expect_that(Fscore(Lw, res$Lw, 1e-1) > .9, is_true())
})


test_that("learn_laplacian_matrix with diamond graph", {
  w <- c(1, 1, 0, 1, 1, 1)
  Lw <- L(w)
  n <- ncol(Lw)
  Y <- MASS::mvrnorm(1000, rep(0, n), MASS::ginv(Lw))
  res <- learn_laplacian_matrix(cov(Y), beta = 40)
  expect_that(relativeError(Lw, res$Lw) < 1e-1, is_true())
  expect_that(Fscore(Lw, res$Lw, 1e-1) > .9, is_true())
})


test_that("learn_bipartite_graph converges with simple bipartite graph", {
  w <- c(1, 0, 0, 1, 0, 1)
  Aw <- A(w)
  n <- ncol(Aw)
  Y <- MASS::mvrnorm(1000, rep(0, n), MASS::ginv(L(w)))
  res <- learn_bipartite_graph(cov(Y))
  expect_that(relativeError(Aw, res$Aw) < 1e-1, is_true())
  expect_that(Fscore(Aw, res$Aw, 1e-1) > .9, is_true())
})


test_that("learn_adjacency_and_laplacian_graph converges with simple bipartite graph", {
  w <- c(1, 0, 0, 1, 0, 1)
  Aw <- A(w)
  n <- ncol(Aw)
  Y <- MASS::mvrnorm(1000, rep(0, n), MASS::ginv(L(w)))
  res <- learn_adjacency_and_laplacian(cov(Y), beta1 = 40)
  expect_that(relativeError(Aw, res$Aw) < 1e-1, is_true())
  expect_that(Fscore(Aw, res$Aw, 1e-1) > .9, is_true())
})


test_that("learn_laplacian_matrix with two components", {
  # test on toy graph from section 3 of https://arxiv.org/pdf/1206.5726.pdf
  Lw <- rbind(c(1, 0, -1, 0),
              c(0, 1, 0, -1),
              c(-1, 0, 1, 0),
              c(0, -1, 0, 1))
  n <- ncol(Lw)
  Y <- MASS::mvrnorm(1000, rep(0, n), MASS::ginv(Lw))
  res <- learn_laplacian_matrix(cov(Y), k = 2, beta = 40)
  expect_that(relativeError(Lw, res$Lw) < 1e-1, is_true())
  expect_that(Fscore(Lw, res$Lw, 1e-1) > .9, is_true())
})


test_that("learn_adjacency_and_laplacian with two components graph", {
  # test on toy graph from section 3 of https://arxiv.org/pdf/1206.5726.pdf
  Lw <- rbind(c(1, 0, -1, 0),
              c(0, 1, 0, -1),
              c(-1, 0, 1, 0),
              c(0, -1, 0, 1))
  n <- ncol(Lw)
  Y <- MASS::mvrnorm(1000, rep(0, n), MASS::ginv(Lw))
  res <- learn_adjacency_and_laplacian(cov(Y), k = 2, beta1 = 20)
  expect_that(relativeError(Lw, res$Lw) < 1e-1, is_true())
  expect_that(Fscore(Lw, res$Lw, 1e-1) > .9, is_true())
})


test_that("learn_laplacian_matrix with two components graph #2", {
  Lw1 <- L(runif(3))
  Lw2 <- L(runif(6))
  n1 <- ncol(Lw1)
  n2 <- ncol(Lw2)

  Lw <- blockDiag(Lw1, Lw2)
  Y <- MASS::mvrnorm(5000, rep(0, n1 + n2), MASS::ginv(Lw))
  res <- learn_laplacian_matrix(cov(Y), k = 2, beta = 40)
  expect_that(relativeError(Lw, res$Lw) < 1e-1, is_true())
  expect_that(Fscore(Lw, res$Lw, 1e-1) > .9, is_true())
})


test_that("learn_adjacency_and_laplacian with two components graph #2", {
  Lw1 <- L(runif(3))
  Lw2 <- L(runif(6))
  n1 <- ncol(Lw1)
  n2 <- ncol(Lw2)

  Lw <- blockDiag(Lw1, Lw2)
  Y <- MASS::mvrnorm(5000, rep(0, n1 + n2), MASS::ginv(Lw))
  res <- learn_adjacency_and_laplacian(cov(Y), k = 2, beta1 = 40, beta2 = 0)
  expect_that(relativeError(Lw, res$Lw) < 1e-1, is_true())
  expect_that(Fscore(Lw, res$Lw, 1e-1) > .9, is_true())
})
