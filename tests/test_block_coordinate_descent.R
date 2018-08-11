library(testthat)
library(spectralGraphTopology)

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
