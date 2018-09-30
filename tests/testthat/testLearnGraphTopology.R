context("testLearnGraphTopology.R")
library(testthat)
library(spectralGraphTopology)

test_that("test_learnGraphTopology_K=1", {
  # test the learning of a single-component graph
  T <- 10000
  w <- sample(1:10, 6)
  K <- 1
  Lw <- L(w)
  N <- ncol(Lw)
  Y <- MASS::mvrnorm(T, as.vector(array(0, N)), MASS::ginv(Lw))
  res <- learnGraphTopology(cov(Y), K, ub=100, beta=.1, maxiter=500)
  expect_that(norm(Lw - res$Lw, type="F") /
              max(1., norm(Lw, type="F")) < 1e-1, is_true())

  # test the learning of a single-component diamond graph
  w <- c(1, 1, 0, 1, 1, 1)
  Lw <- L(w)
  Y <- MASS::mvrnorm(T, as.vector(array(0, N)), MASS::ginv(Lw))
  res <- learnGraphTopology(cov(Y), K, ub=10, beta=10, maxiter=500)
  expect_that(norm(Lw - res$Lw, type="F") /
              max(1., norm(Lw, type="F")) < 1e-1, is_true())
})


# test on toy graph from section 3 of https://arxiv.org/pdf/1206.5726.pdf
test_that("test_learnGraphTopology_K=2", {
  T <- 10000
  K <- 2
  N <- 4
  Lw <- rbind(c(1, 0, -1, 0),
              c(0, 1, 0, -1),
              c(-1, 0, 1, 0),
              c(0, -1, 0, 1))
  Y <- MASS::mvrnorm(T, as.vector(array(0, N)), MASS::ginv(Lw))
  res <- learnGraphTopology(cov(Y), K, beta=20.)
  expect_that(norm(Lw - res$Lw, type="F") / norm(Lw, type="F") < 1e-1,
              is_true())
})


test_that("test_learnGraphTopology_K=2", {
  T <- 2000
  w1 <- runif(3)
  w2 <- runif(6)
  K <- 2
  Lw1 <- L(w1)
  Lw2 <- L(w2)
  N1 <- ncol(Lw1)
  N2 <- ncol(Lw2)

  Lw <- blockDiag(Lw1, Lw2)
  Y <- MASS::mvrnorm(T, rep(0, N1 + N2), MASS::ginv(Lw))
  res <- learnGraphTopology(cov(Y), K, beta=100)
  expect_that(norm(Lw - res$Lw, type="F") / norm(Lw, type="F") < 1e-1,
              is_true())
})
