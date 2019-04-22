context("Laplacian and Adjacency matrices estimation")
library(testthat)
library(spectralGraphTopology)


test_that("learn_laplacian_matrix with single component random graph", {
  w <- sample(1:10, 6)
  Laplacian <- L(w)
  n <- ncol(Laplacian)
  Y <- MASS::mvrnorm(n * 500, rep(0, n), MASS::ginv(Laplacian))
  res <- learn_laplacian_matrix(cov(Y))
  expect_that(res$convergence, is_true())
  expect_that(relative_error(Laplacian, res$Laplacian) < 1e-1, is_true())
  expect_that(metrics(Laplacian, res$Laplacian, 1e-1)[1] > .9, is_true())
})


test_that("learn_laplacian_matrix with diamond graph", {
  w <- c(1, 1, 0, 1, 1, 1)
  Laplacian <- L(w)
  n <- ncol(Laplacian)
  Y <- MASS::mvrnorm(n * 500, rep(0, n), MASS::ginv(Laplacian))
  res <- learn_laplacian_matrix(cov(Y))
  expect_that(res$convergence, is_true())
  expect_that(relative_error(Laplacian, res$Laplacian) < 1e-1, is_true())
  expect_that(metrics(Laplacian, res$Laplacian, 1e-1)[1] > .9, is_true())
})


test_that("learn_bipartite_graph converges with simple bipartite graph", {
  w <- c(1, 0, 0, 1, 0, 1)
  Adjacency <- A(w)
  n <- ncol(Adjacency)
  Y <- MASS::mvrnorm(n * 500, rep(0, n), MASS::ginv(L(w)))
  res <- learn_bipartite_graph(cov(Y))
  expect_that(res$convergence, is_true())
  expect_that(relative_error(Adjacency, res$Adjacency) < 1e-1, is_true())
  expect_that(metrics(Adjacency, res$Adjacency, 1e-1)[1] > .9, is_true())
})


test_that("learn_adjancecy_and_laplacian can learn k-component bipartite graph", {
  w <- c(1, 0, 0, 1, 0, 1)
  Laplacian <- block_diag(L(w), L(w))
  n <- ncol(Laplacian)
  Y <- MASS::mvrnorm(2 * n * 200, rep(0, n), MASS::ginv(Laplacian))
  graph <- learn_adjacency_and_laplacian(cov(Y), k = 2)
  expect_that(graph$convergence, is_true())
  expect_that(relative_error(Laplacian, graph$Laplacian) < 1e-1, is_true())
  expect_that(metrics(Laplacian, graph$Laplacian, 1e-1)[1] > .9, is_true())
})


test_that("learn_adjacency_and_laplacian_graph converges with simple bipartite graph", {
  w <- c(1, 0, 0, 1, 0, 1)
  Adjacency <- A(w)
  n <- ncol(Adjacency)
  Y <- MASS::mvrnorm(n * 500, rep(0, n), MASS::ginv(L(w)))
  res <- learn_adjacency_and_laplacian(cov(Y))
  expect_that(res$convergence, is_true())
  expect_that(relative_error(Adjacency, res$Adjacency) < 1e-1, is_true())
  expect_that(metrics(Adjacency, res$Adjacency, 1e-1)[1] > .9, is_true())
})


test_that("learn_laplacian_matrix with two components", {
  Laplacian <- rbind(c(1, -1, 0, 0),
                     c(-1, 1, 0, 0),
                     c(0, 0, 1, -1),
                     c(0, 0, -1, 1))
  n <- ncol(Laplacian)
  Y <- MASS::mvrnorm(n * 500, rep(0, n), MASS::ginv(Laplacian))
  res <- learn_laplacian_matrix(cov(Y), k = 2, beta = 10)
  expect_that(res$convergence, is_true())
  expect_that(relative_error(Laplacian, res$Laplacian) < 1e-1, is_true())
  expect_that(metrics(Laplacian, res$Laplacian, 1e-1)[1] > .9, is_true())
})


test_that("check that learn_adjacency_and_laplacian and learn_laplacian_matrix are consistent", {
  w <- c(1, 0, 0, 1, 0, 1)
  Laplacian <- L(w)
  n <- ncol(Laplacian)
  Y <- MASS::mvrnorm(n * 100, rep(0, n), MASS::ginv(Laplacian))
  S <- cov(Y)
  res1 <- learn_adjacency_and_laplacian(S, w0 = "qp", beta = 10, nu = 0)
  res2 <- learn_laplacian_matrix(S, w0 = "qp", beta = 10)
  expect_that(res1$convergence, is_true())
  expect_that(res2$convergence, is_true())
  expect_that(all(abs(res1$obj_fun - res2$obj_fun) < 1e-9), is_true())
  expect_that(all(abs(res1$w - res2$w) < 1e-9), is_true())
})


test_that("learn_laplacian_matrix with two components graph #2", {
  Laplacian1 <- L(runif(3))
  Laplacian2 <- L(runif(6))
  n1 <- ncol(Laplacian1)
  n2 <- ncol(Laplacian2)

  Laplacian <- block_diag(Laplacian1, Laplacian2)
  Y <- MASS::mvrnorm(500 * (n1 + n2), rep(0, n1 + n2), MASS::ginv(Laplacian))
  res <- learn_laplacian_matrix(cov(Y), k = 2)
  expect_that(res$convergence, is_true())
  expect_that(relative_error(Laplacian, res$Laplacian) < 1e-1, is_true())
  expect_that(metrics(Laplacian, res$Laplacian, 1e-1)[1] > .85, is_true())
})


test_that("learn_adjacency_and_laplacian with two components graph #2", {
  Laplacian1 <- L(runif(3))
  Laplacian2 <- L(runif(6))
  n1 <- ncol(Laplacian1)
  n2 <- ncol(Laplacian2)

  Laplacian <- block_diag(Laplacian1, Laplacian2)
  Y <- MASS::mvrnorm(500 * (n1 + n2), rep(0, n1 + n2), MASS::ginv(Laplacian))
  res <- learn_adjacency_and_laplacian(cov(Y), k = 2, w0 = "qp", nu = 0, eig_tol = 1e-5)
  expect_that(res$convergence, is_true())
  expect_that(relative_error(Laplacian, res$Laplacian) < 1e-1, is_true())
  expect_that(metrics(Laplacian, res$Laplacian, 1e-1)[1] > .9, is_true())
})


test_that("learn_dregular_graph works as expected", {
   Lw1 <- L(c(1, 1, 1))
   Lw2 <- L(c(1, 0, 1, 1, 0, 1))
   Lw <- block_diag(Lw1, Lw2)
   n1 <- ncol(Lw1)
   n2 <- ncol(Lw2)
   Y <- MASS::mvrnorm(500 * (n1 + n2), rep(0, n1 + n2), MASS::ginv(Lw))
   res <- learn_dregular_graph(cov(Y), k = 2, w0 = "qp")
   expect_that(res$convergence, is_true())
   expect_that(relative_error(Lw, res$Laplacian) < 1e-1, is_true())
   expect_that(metrics(Lw, res$Laplacian, 1e-1)[1] > .9, is_true())
})

with_parameters_test_that("learn_normalized_laplacian works on toy data", {
  n <- ncol(Laplacian)
  Y <- MASS::mvrnorm(n * 5, rep(0, n), MASS::ginv(Laplacian))
  res <- learn_normalized_laplacian(cov(Y), scale = FALSE, k = 2, record_objective = TRUE, reltol = 1e-3)
  expect_that(res$convergence, is_true())
  expect_that(relative_error(Laplacian, res$NormalizedLaplacian) < 1e-1, is_true())
  expect_that(metrics(Laplacian, res$NormalizedLaplacian, 1e-1)[1] > .9, is_true())
}, cases(list(Laplacian = rbind(c( 1, -1,  0,  0),
                                c(-1,  1,  0,  0),
                                c( 0,  0,  1, -1),
                                c( 0,  0, -1,  1))),
         list(Laplacian = rbind(c(   1, -1/2, -1/2,  0,  0),
                                c(-1/2,    1, -1/2,  0,  0),
                                c(-1/2, -1/2,    1,  0,  0),
                                c(   0,    0,    0,  1, -1),
                                c(   0,    0,    0, -1,  1)))
        )
)

test_that("learn_normalized_laplacian works on non-normalized data", {
  Laplacian1 <- L(runif(3))
  Laplacian2 <- L(runif(6))
  n1 <- ncol(Laplacian1)
  n2 <- ncol(Laplacian2)

  Laplacian <- block_diag(Laplacian1, Laplacian2)
  Y <- MASS::mvrnorm(5 * (n1 + n2), rep(0, n1 + n2), MASS::ginv(Laplacian))
  res <- learn_normalized_laplacian(cov(Y), k = 2)
  expect_that(res$convergence, is_true())
  expect_that(metrics(Laplacian, res$NormalizedLaplacian, 1e-1)[1] > .9, is_true())
})
