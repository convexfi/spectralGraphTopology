library(testthat)
library(spectralGraphTopology)

# lambda update step using CVX for the
# sake of unit testing
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
  U <- U_update(w, N, K)
  expect_that(all.equal(crossprod(U), diag(array(1., N-K)),
                         check.attributes = FALSE), is_true())
  expect_that(ncol(U) == N-K, is_true())
  expect_that(nrow(U) == N, is_true())
})

test_that("test_lambda_update_consistency", {
  # test that the eigen values meet the criterion after being updated
  w <- runif(1000)
  N <- as.integer(.5 * (1 + sqrt(1 + 8 * length(w))))
  K <- 1
  U <- U_update(w, N, K)
  lb <- 1e-2
  ub <- 100.
  beta <- .5
  q <- N - K

  lambda <- lambda_update(lb, ub, beta, U, w, N, K)
  expect_that(length(lambda) == q, is_true())
  expect_that(all(lambda[1] >= lb, lambda[q] <= ub,
                  lambda[2:q] >= lambda[1:(q-1)]), is_true())

  # compare against results from CVXR
  lambda_cvx <- lambda_update_cvx(lb, ub, beta, U, w, N, K)
  expect_that(length(lambda_cvx) == q, is_true())
  expect_that(all(lambda_cvx[1] >= lb, lambda_cvx[q] <= ub,
                  lambda_cvx[2:q] >= lambda_cvx[1:(q-1)]), is_true())

  expect_that(all(abs(lambda_cvx - lambda) < 1e-3), is_true())
})

test_that("test_lambda_update_equals_to_ub", {
  w <- runif(6)
  N <- as.integer(.5 * (1 + sqrt(1 + 8 * length(w))))
  K <- 1
  U <- U_update(w, N, K)
  lb <- 1e-2
  ub <- 1.5
  beta <- .5
  q <- N - K

  lambda <- lambda_update(lb, ub, beta, U, w, N, K)
  lambda_cvx <- lambda_update_cvx(lb, ub, beta, U, w, N, K)
  expect_that(all(abs(lambda_cvx - lambda) < 1e-3), is_true())

  expect_that(length(lambda) == q, is_true())
  expect_that(all(lambda[1] >= lb, lambda[q] <= ub,
                  lambda[2:q] >= lambda[1:(q-1)]), is_true())

  expect_that(length(lambda_cvx) == q, is_true())
  expect_that(all(lambda_cvx[1] >= lb, lambda_cvx[q] <= ub,
                  lambda_cvx[2:q] >= lambda_cvx[1:(q-1)]), is_true())
})
