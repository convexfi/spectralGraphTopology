library(testthat)
library(spectralGraphTopology)

test_that("test_learn_graph_topology", {
  n <- 4
  k <- 10000
  w <- c(1:6)
  K <- 1
  Lw <- LOp(w, n)
  Y <- t(MASS::mvrnorm(k, as.vector(array(0, n)), MASS::ginv(Lw)))
  Lw_est <- learnGraphTopology(Y, K, alpha=0, ub=1e4, beta=.1, pho=0, maxiter=200)
  expect_that(norm(Lw - Lw_est, type="F") / max(1., norm(Lw, type="F")) < 1e-1, is_true())
})
