# this module contains implementations of graph learning methods
# from smooth signals

#' @title Learn graphs from a smooth signal representation approach
#'
#' This function learns a graph from a observed data matrix using the
#' method proposed by Dong (2016).
#'
#' @references X. Dong, D. Thanou, P. Frossard and P. Vandergheynst, "Learning
#'             Laplacian Matrix in Smooth Graph Signal Representations,"
#'             in IEEE Transactions on Signal Processing, vol. 64, no. 23,
#'             pp. 6160-6173, Dec.1, 2016.
#' @export
learn_graph_sigrep <- function(X, alpha = 1e-2, beta = 1e-1, maxiter = 1000, ftol = 1e-4) {
  n <- nrow(X)
  Y <- X
  obj_values <- c()
  fun0 <- Inf
  if (verbose)
    pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta  relerr: :relerr",
                                     total = maxiter, clear = FALSE, width = 80)
  for (k in c(1:maxiter)) {
    L <- glsigrep.update_L(Y, alpha, beta, n)
    Y <- glsigrep.update_Y(L, X, alpha, n)
    funk <- glsigrep.obj_function(X, Y, L, alpha, beta)
    obj_values <- c(obj_values, funk)
    relerr <- abs(funk - fun0) / fun0
    if (k > 1 && abs(funk - fun0) / fun0 < ftol) break
    if (verbose) pb$tick(token = list(relerr = relerr))
    fun0 <- funk
  }
  return(list(L, Y, obj_values))
}


glsigrep.update_L <- function(Y, alpha, beta, n) {
  L <- CVXR::Symmetric(n, n)
  ones <- rep(1, n)
  zeros <- rep(0, n)
  zeros_mat <- matrix(0, n, n)
  objective <- CVXR::Minimize(alpha * CVXR::matrix_trace(L %*% Y %*% t(Y))
                              + beta * CVXR::sum_squares(L))
  constraints <- list(CVXR::matrix_trace(L) == n,
                      L == t(L),
                      L %*% ones == zeros,
                      L - CVXR::diag(CVXR::diag(L)) <= zeros_mat)
  problem <- CVXR::Problem(objective, constraints)
  result <- solve(problem)
  return(result$getValue(L))
}

glsigrep.update_Y <- function(L, X, alpha, n) {
  In <- diag(n)
  return(solve(In + alpha * L, X))
}

glsigrep.obj_function <- function(X, Y, L, alpha, beta) {
  return(norm(X - Y, 'F') ^ 2 + alpha * sum(diag(L %*% Y %*% t(Y))) + beta * norm(L, 'F') ^ 2)
}
