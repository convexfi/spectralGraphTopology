library(testthat)
library(spectralGraphTopology)

test_that("test_learnGraphTopologyK=1", {
  # test the learning of a single-component graph
  T <- 10000
  w <- sample(1:10, 6)
  N <- as.integer(.5 * (1 + sqrt(1 + 8 * length(w))))
  K <- 1
  Theta <- L(w)
  Y <- t(MASS::mvrnorm(T, as.vector(array(0, N)), MASS::ginv(Theta)))
  Theta_est <- learnGraphTopology(Y, K, ub=100, beta=.1, maxiter=500)
  expect_that(norm(Theta - Theta_est, type="F") /
              max(1., norm(Theta, type="F")) < 1e-1, is_true())

  # test the learning of a single-component diamond graph
  w <- c(1, 1, 0, 1, 1, 1)
  Theta <- L(w)
  Y <- t(MASS::mvrnorm(T, as.vector(array(0, N)), MASS::ginv(Theta)))
  Theta_est <- learnGraphTopology(Y, K, ub=10, beta=10, maxiter=500)
  expect_that(norm(Theta - Theta_est, type="F") /
              max(1., norm(Theta, type="F")) < 1e-1, is_true())
})


# This test is currently failing, it should pass after we make K > 1 work
test_that("test_learnGraphTopology_K=2", {
  T <- 10000
  w1 <- runif(3)
  w2 <- runif(6)
  K <- 2
  Theta1 <- L(w1)
  Theta2 <- L(w2)
  N1 <- as.integer(.5 * (1 + sqrt(1 + 8 * length(w1))))
  N2 <- as.integer(.5 * (1 + sqrt(1 + 8 * length(w2))))

  Theta <- rbind(cbind(Theta1, matrix(0, N1, N2)),
                 cbind(matrix(0, N2, N1), Theta2))
  Y <- t(MASS::mvrnorm(T, as.vector(array(0, N1 + N2)), MASS::ginv(Theta)))
  Theta_est <- learnGraphTopology(Y, K, ub=100, beta=1., maxiter=500)
  expect_that(norm(Theta - Theta_est, type="F") /
              max(1., norm(Theta, type="F")) < 1e-1, is_true())
})
