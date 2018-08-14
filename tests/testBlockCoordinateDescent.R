library(testthat)
library(spectralGraphTopology)


test_that("test_U_update_consistency", {
  # test that U remains orthonormal after being updated
  w <- c(1, 2, 3, 4, 5, 6)
  N <- 4
  K <- 1
  U <- U_update(w, N, K)
  expect_that(all.equal(crossprod(U), diag(array(1., N-K)),
                         check.attributes = FALSE), is_true())
  expect_that(ncol(U) == N-K, is_true())
  expect_that(nrow(U) == N, is_true())
})

test_that("test_lambda_update_consistency", {
  # test that the eigen values meet the criterion after being updated
  w <- runif(6)
  N <- 4
  K <- 1
  U <- U_update(w, N, K)
  lb <- 1e-2
  ub <- 100.
  beta <- .5
  l <- N - K

  lambda <- lambda_update(lb, ub, beta, U, w, N, K)
  expect_that(all(lambda[l] >= lb, lambda[1] <= ub,
                  lambda[1:(l-1)] >= lambda[2:l]), is_true())
})
