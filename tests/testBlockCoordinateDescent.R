library(testthat)
library(spectralGraphTopology)

# lambda update step using CVX for the sake of unit testing
lambda_update_cvx <- function(lb, ub, beta, U, w, N, K) {
  d <- diag(t(U) %*% L(w) %*% U)
  q <- N - K
  lambda <- CVXR::Variable(q)
  objective <- CVXR::Minimize(sum(.5 * beta * (lambda - d)^2 - log(lambda)))
  constraints <- list(lambda[q] <= ub, lambda[1] >= lb, lambda[2:q] >= lambda[1:(q-1)])
  prob <- CVXR::Problem(objective, constraints)
  result <- solve(prob)
  return(as.vector(result$getValue(lambda)))
}

test_that("test_U_update_consistency", {
  # test that U remains orthonormal after being updated
  w <- runif(1000)
  N <- as.integer(.5 * (1 + sqrt(1 + 8 * length(w))))
  K <- 1
  U <- U_update(w, N)
  expect_that(all.equal(crossprod(U), diag(array(1., N)),
                         check.attributes = FALSE), is_true())
  expect_that(ncol(U) == N, is_true())
  expect_that(nrow(U) == N, is_true())
})

test_that("test_lambda_update_consistency", {
  # test that the eigen values meet the criterion after being updated
  w <- runif(1000)
  N <- as.integer(.5 * (1 + sqrt(1 + 8 * length(w))))
  K <- 1
  U <- U_update(w, N)
  lb <- 1e-2
  ub <- 10
  beta <- .5
  q <- N-K

  lambda <- lambda_update(lb, ub, beta, U, w, N, K)
  expect_that(length(lambda) == N, is_true())
  expect_that(lambda[K] == 0, is_true())
  expect_that(all(lambda[K+1] >= lb, lambda[N] <= ub,
                  lambda[(K+2):N] >= lambda[(K+1):(N-1)]), is_true())

  # compare against results from CVXR
  lambda_cvx <- lambda_update_cvx(lb, ub, beta, U[, (K+1):N], w, N, K)
  expect_that(length(lambda_cvx) == q, is_true())
  expect_that(all(abs(lambda_cvx - lambda[(K+1):N]) < 1e-3), is_true())
})

test_that("test_lambda_update_equals_to_ub", {
  w <- runif(6)
  N <- as.integer(.5 * (1 + sqrt(1 + 8 * length(w))))
  K <- 1
  U <- U_update(w, N)
  lb <- 1e-2
  ub <- 1.5
  beta <- .5

  lambda <- lambda_update(lb, ub, beta, U, w, N, K)
  expect_that(length(lambda) == N, is_true())
  expect_that(lambda[K] == 0, is_true())
  expect_that(all(lambda[K+1] >= lb, lambda[N] <= ub,
                  lambda[(K+2):N] >= lambda[(K+1):(N-1)]), is_true())
  lambda_cvx <- lambda_update_cvx(lb, ub, beta, U[, (K+1):N], w, N, K)
  expect_that(all(abs(lambda_cvx - lambda[(K+1):N]) < 1e-4), is_true())
})

test_that("test_lambda_update_equals_to_lb", {
  w <- runif(20)
  N <- as.integer(.5 * (1 + sqrt(1 + 8 * length(w))))
  K <- 1
  U <- U_update(w, N)
  lb <- 3
  ub <- 100
  beta <- .5

  lambda <- lambda_update(lb, ub, beta, U, w, N, K)
  expect_that(length(lambda) == N, is_true())
  expect_that(lambda[K] == 0, is_true())
  expect_that(all(lambda[K+1] >= lb, lambda[N] <= ub,
                  lambda[(K+2):N] >= lambda[(K+1):(N-1)]), is_true())
  lambda_cvx <- lambda_update_cvx(lb, ub, beta, U[, (K+1):N], w, N, K)
  expect_that(all(abs(lambda_cvx - lambda[(K+1):N]) < 1e-4), is_true())
})


test_that("test_lambda_update_with_terms_equal_to_both_lb_and_ub", {
  w <- runif(20)
  N <- as.integer(.5 * (1 + sqrt(1 + 8 * length(w))))
  K <- 1
  U <- U_update(w, N)
  lb <- 3.3
  ub <- 4.
  beta <- .5

  lambda <- lambda_update(lb, ub, beta, U, w, N, K)
  expect_that(length(lambda) == N, is_true())
  expect_that(lambda[K] == 0, is_true())
  expect_that(all(lambda[K+1] >= lb, lambda[N] <= ub,
                  lambda[(K+2):N] >= lambda[(K+1):(N-1)]), is_true())
  lambda_cvx <- lambda_update_cvx(lb, ub, beta, U[, (K+1):N], w, N, K)
  expect_that(all(abs(lambda_cvx - lambda[(K+1):N]) < 1e-4), is_true())
})
