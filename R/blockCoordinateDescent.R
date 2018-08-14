# This module defines the updates of the block coordinate descent algorithm.
# More precisely, there are three variables to be updated:
#  1) a vector w (see expression (9)),
#  2) a matrix U (see expression (11)),
#  3) a matrix Λ (see expression (15)).

w_update <- function(w, U, beta, lambda, N, Km) {
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
  #grad_f <- Lstar(L(w) - tcrossprod(U*sqrt(lambda)) + Km / beta)
  grad_f <- Lstar(L(w) - U %*% diag(lambda) %*% t(U) + Km / beta)
  w_update <- w - .5 * grad_f / N
  return(pmax(0, w_update))
}


U_update <- function(w, N, K) {
  # Updates the value of U
  #
  # Args:
  #   w: vector
  #   N: dimension of each data sample
  #   K: number of components
  #
  # Returns:
  #   U_update: the updated value of U
  return(eigen(L(w))$vectors[, 1:(N-K)])
}


lambda_update <- function(lb, ub, beta, U, w, N, K) {
  # Updates the value of U
  #
  # Args:
  #   lb, ub: lower and upper bounds on the Lambda vector components
  #   N: dimension of each data sample
  #   K: number of components of the graph
  #
  # Returns:
  #   Lambda_update: the updated value of Lambda

  d <- diag(t(U) %*% L(w) %*% U)
  l <- N - K
  lambda <- CVXR::Variable(l)
  objective <- CVXR::Minimize(sum(.5 * beta * (lambda - d)^2 - log(lambda)))
  constraints <- list(lambda[l] >= lb, lambda[1] <= ub,
                      lambda[1:(l-1)] >= lambda[2:l])
  prob <- CVXR::Problem(objective, constraints)
  prob_data <- CVXR::get_problem_data(prob, solver = "ECOS")
  solver_output <- ECOSolveR::ECOS_csolve(c = prob_data[["c"]],
                                            G = prob_data[["G"]],
                                            h = prob_data[["h"]],
                                            dims = prob_data[["dims"]],
                                            A = prob_data[["A"]],
                                            b = prob_data[["b"]])
  result <- CVXR::unpack_results(prob, "ECOS", solver_output)
  return(as.vector(result$getValue(lambda)))
}
