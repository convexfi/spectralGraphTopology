context("testBlockCoordinateDescent.R")
library(testthat)
library(patrick)
library(spectralGraphTopology)

# lambda update step using CVX for the sake of unit testing
lambda_update_cvx <- function(lb, ub, beta, U, Lw, N, K) {
  d <- diag(t(U) %*% Lw %*% U)
  q <- N - K
  lambda <- CVXR::Variable(q)
  objective <- CVXR::Minimize(sum(.5 * beta * (lambda - d)^2 - log(lambda)))
  constraints <- list(lambda[q] <= ub, lambda[1] >= lb, lambda[2:q] >= lambda[1:(q-1)])
  prob <- CVXR::Problem(objective, constraints)
  result <- solve(prob)
  return(as.vector(result$getValue(lambda)))
}


psi_update_cvx <- function(V, Aw) {
  c <- diag(t(V) %*% Aw %*% V)
  q <- length(c)
  psi <- CVXR::Variable(q)
  objective <- CVXR::Minimize(sum((psi - c) ^ 2))
  constraints <- list(psi[1:(q - 1)] < psi[2:q])
  j <- length(constraints) + 1
  i <- 1
  while (i < q + 1 - i) {
    constraints[[j]] <- psi[i] == -psi[q + 1 - i]
    j <- j + 1
    i <- i + 1
  }
  prob <- CVXR::Problem(objective, constraints)
  result <- solve(prob)
  return(as.vector(result$getValue(psi)))
}


test_that("test that U remains orthonormal after being updated", {
  w <- runif(1000)
  N <- as.integer(.5 * (1 + sqrt(1 + 8 * length(w))))
  K <- 1
  U <- laplacian.U_update(L(w), N, K)
  q <- N - K
  expect_that(all.equal(crossprod(U), diag(array(1., q)),
                        check.attributes = FALSE), is_true())
  expect_that(ncol(U) == q, is_true())
  expect_that(nrow(U) == N, is_true())
})


test_that("test that V remains orthonormal after being updated", {
  w <- runif(4*9)
  n <- as.integer(.5 * (1 + sqrt(1 + 8 * length(w))))
  z <- 3
  V <- adjacency.V_update(A(w), n, z)
  q <- n - z
  expect_that(all.equal(crossprod(V), diag(array(1., q)),
                        check.attributes = FALSE), is_true())
  expect_that(ncol(V) == q, is_true())
  expect_that(nrow(V) == n, is_true())
})


test_that("test that the eigenvalues of the adjacency matrix meet the criterion", {
  w <- runif(4*9)
  n <- as.integer(.5 * (1 + sqrt(1 + 8 * length(w))))
  z <- 3
  Aw <- A(w)
  V <- adjacency.V_update(Aw, n, z)
  psi <- adjacency.psi_update(V, Aw)
  psi_cvx <- psi_update_cvx(V, Aw)
  expect_equal(psi, psi_cvx, tolerance = 1e-4)
})


with_parameters_test_that("test that the eigenvalues of the Laplacian matrix
                          meet the criterion after being updated", {
    N <- as.integer(.5 * (1 + sqrt(1 + 8 * length(w))))
    K <- 1
    Lw <- L(w)
    U <- laplacian.U_update(Lw, N, K)
    beta <- .5
    q <- N - K

    lambda <- laplacian.lambda_update(lb, ub, beta, U, Lw, N, K)
    expect_that(length(lambda) == q, is_true())
    expect_that(all(lambda[1] >= lb, lambda[q] <= ub,
                    lambda[2:q] >= lambda[1:(q-1)]), is_true())

    # compare against results from CVXR
    lambda_cvx <- lambda_update_cvx(lb, ub, beta, U, Lw, N, K)
    expect_that(length(lambda_cvx) == q, is_true())
    expect_that(all(abs(lambda_cvx - lambda) < 1e-3), is_true())
  },
  cases(
        list(lb = 1e-2, ub = 10,  w = runif(1000)),
        list(lb = 1e-2, ub = 1.5, w = runif(6)),
        list(lb = 3,    ub = 100, w = runif(20)),
        list(lb = 3.3,  ub = 4,   w = runif(20))
       )
)
