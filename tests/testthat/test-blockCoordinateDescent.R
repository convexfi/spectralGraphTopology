context("block coordinate descent")
library(testthat)
library(patrick)
library(spectralGraphTopology)

# lambda update step using CVX for the sake of unit testing
lambda_update_cvx <- function(lb, ub, beta, U, Lw, k) {
  n <- ncol(Lw)
  d <- diag(t(U) %*% Lw %*% U)
  q <- n - k
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
  n <- as.integer(.5 * (1 + sqrt(1 + 8 * length(w))))
  k <- 1
  U <- laplacian.U_update(L(w), k)
  q <- n - k
  expect_true(all.equal(crossprod(U), diag(array(1., q)),
                        check.attributes = FALSE))
  expect_true(ncol(U) == q)
  expect_true(nrow(U) == n)
})


test_that("test that V remains orthonormal after being updated", {
  w <- runif(4*9)
  n <- as.integer(.5 * (1 + sqrt(1 + 8 * length(w))))
  z <- 3
  V <- bipartite.V_update(A(w), z)
  q <- n - z
  expect_true(all.equal(crossprod(V), diag(array(1., q)),
                        check.attributes = FALSE))
  expect_true(ncol(V) == q)
  expect_true(nrow(V) == n)
})


test_that("test that the eigenvalues of the adjacency matrix meet the criterion", {
  skip_if_not_installed("CVXR")
  w <- runif(4*9)
  n <- as.integer(.5 * (1 + sqrt(1 + 8 * length(w))))
  z <- 3
  Aw <- A(w)
  V <- bipartite.V_update(Aw, z)
  psi <- bipartite.psi_update(V, Aw)
  psi_cvx <- psi_update_cvx(V, Aw)
  expect_equal(psi, psi_cvx, tolerance = 1e-4)
})


with_parameters_test_that("test that the eigenvalues of the Laplacian matrix
                          meet the criterion after being updated", {
    skip_if_not_installed("CVXR")
    n <- as.integer(.5 * (1 + sqrt(1 + 8 * length(w))))
    k <- 1
    Lw <- L(w)
    U <- laplacian.U_update(Lw, k)
    beta <- .5
    q <- n - k

    lambda <- laplacian.lambda_update(lb, ub, beta, U, Lw, k)
    expect_true(length(lambda) == q)
    expect_true(all(lambda[1] >= lb, lambda[q] <= ub,
                    lambda[2:q] >= lambda[1:(q-1)]))

    # compare against results from CVXR
    lambda_cvx <- lambda_update_cvx(lb, ub, beta, U, Lw, k)
    expect_true(length(lambda_cvx) == q)
    expect_true(all(abs(lambda_cvx - lambda) < 1e-3))
  },
  cases(
        list(lb = 1e-2, ub = 10,  w = runif(1000)),
        list(lb = 1e-2, ub = 1.5, w = runif(6)),
        list(lb = 3,    ub = 100, w = runif(20)),
        list(lb = 3.3,  ub = 4,   w = runif(20))
       )
)
