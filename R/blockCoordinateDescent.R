# This module defines the updates of the block coordinate descent algorithm.
# More precisely, there are three variables to be updated:
#  1) a vector w (see expression (9)),
#  2) a matrix U (see expression (11)),
#  3) a matrix Λ (see expression (15)).

w_update <- function(w, U, beta, Lambda, n, Km) {
  # Updates the value of w
  #
  # Args:
  #   w: vector
  #   U: matrix
  #   β: scalar
  #   Λ: vector
  #   n: dimension of each data sample
  #   Km: matrix (see section 1.1 for its definition)
  #
  # Returns:
  #   w_update: the updated value of w
  grad_f <- LStarOp(LOp(w, n)) - LStarOp(U %*% diag(Lambda) %*% t(U) - Km / beta)
  w_update <- w - .5 * grad_f / n
  mask <- w_update < 0
  w_update[mask] <- 0
  return(w_update)
}


U_update <- function(w, n, K) {
  # Updates the value of U
  #
  # Args:
  #   w: vector
  #   n: dimension of each data sample
  #   K: number of components
  #
  # Returns:
  #   U_update: the updated value of U
  return(eigen(LOp(w, n))$vectors[, 1:(n-K)])
}


Lambda_update <- function(lb, ub, beta, U, w, n, K) {
  # Updates the value of U
  #
  # Args:
  #   lb, ub: lower and upper bounds on the Lambda vector components
  #   n: dimension of each data sample
  #   K: number of components of the graph
  #
  # Returns:
  #   Lambda_update: the updated value of Lambda

  d <- diag(t(U) %*% LOp(w, n) %*% U)
  l <- n - K
  lambda <- CVXR::Variable(l)
  objective <- CVXR::Minimize(sum(.5 * beta * (lambda - d)^2 - log(lambda)))
  constraints <- list(lambda[l] >= lb, lambda[1] <= ub,
                      lambda[1:(l-1)] >= lambda[2:l])
  prob <- CVXR::Problem(objective, constraints)
  result <- solve(prob)
  return(as.vector(result$getValue(lambda)))
}
