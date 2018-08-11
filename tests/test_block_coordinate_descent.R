library(testthat)
library(spectralGraphTopology)

test_that("test_w_update_runs", {
  n <- 4
  w <- c(1, 2, 3, 4, 5, 6)
  K <- 1
  Lw <- LOp(w, n)
  Y <- MASS::mvrnorm(n, c(0, 0, 0, 0), MASS::ginv(Lw))
  S <- Y %*% t(Y) / 4
  alpha <- 1.
  H <- alpha * (2. * diag(n) - matrix(1, n, n))
  Km <- S + H
  lb <- 1e-2
  ub <- 10.
  beta <- .5
  w_ <- w+.1*runif(6)
  U <- U_update(w_, n, K)
  Lambda <- Lambda_update(lb, ub, beta, U, w_, n, K)
  w <- w_update(w_, U, beta, Lambda, n, Km)
})

test_that("test_U_update_consistency", {
  w <- c(1, 2, 3, 4, 5, 6)
  n <- 4
  K <- 1
  U <- U_update(w, n, K)
  # test that U is orthonormal
  expect_that(all.equal(t(U) %*% U, diag(array(1., n-K)),
                         check.attributes = FALSE), is_true())
})

test_that("test_Lambda_update_consistency", {
  w <- runif(6)
  n <- 4
  K <- 1
  U <- U_update(w, n, K)
  lb <- 1e-2
  ub <- 2.
  beta <- .5
  l <- n - K

  lambda <- Lambda_update(lb, ub, beta, U, w, n, K)
  expect_that(all(lambda[1] >= lb, lambda[l] <= ub,
                  lambda[1:(l-1)] <= lambda[2:l]), is_true())
})
