library(testthat)
library(spectralGraphTopology)

test_that("test_learn_graph_topology", {
  n <- 4
  k <- 1000
  w <- runif(as.integer(.5 * n * (n - 1)))
  K <- 1
  Lw <- LOp(w, n)
  Y <- t(MASS::mvrnorm(k, as.vector(array(0, n)), MASS::ginv(Lw)))
  Lw_est <- learnGraphTopology(Y, K)
  expect_that(norm(Lw - Lw_est, type="F") / max(1., norm(Lw, type="F")) < 1e-4, is_true())
})
