# this module contains implementations of the algorithms
# proposed in L. Zhao et. al. "Optimization Algorithms for
# Graph Laplacian Estimation via ADMM and MM", IEEE Trans. Sign. Proc. 2019.

get_incidence_from_adjacency <- function(A) {
  p <- nrow(A)
  #non_zero_edges <- sum(Ainv(A) > 0)
  E <- matrix(0, p, .5 * p * (p-1))
  k <- 1
  for (i in c(1:(p-1))) {
    for (j in c((i+1):p)) {
      if (A[i, j] > 0) {
        E[i, k] <- 1
        E[j, k] <- -1
        k <- k + 1
      }
    }
  }
  return(E)
}

#' @title Learn the weighted Laplacian matrix of a graph given the topological connections
#'
#' @param S a pxp sample covariance/correlation matrix
#' @param A the binary adjacency matrix of the graph
#' @param w0 initial estimate for the weight vector the graph or a string
#'        selecting an appropriate method. Available methods are: "qp": finds w0 that minimizes
#'        ||ginv(S) - L(w0)||_F, w0 >= 0; "naive": takes w0 as the negative of the
#'        off-diagonal elements of the pseudo inverse, setting to 0 any elements s.t.
#'        w0 < 0
#' @param alpha L1 regularization hyperparameter
#' @param maxiter the maximum number of iterations
#' @param abstol absolute tolerance on the weight vector w
#' @param reltol relative tolerance on the weight vector w

learn_laplacian_gle_mm <- function(S, A, w0 = "naive", alpha = 0, maxiter = 1000,
                                   reltol = 1e-4, abstol = 1e-5, record_objective = FALSE) {
  Sinv <- MASS::ginv(S)
  w <- w_init(w0, Sinv)
  wk <- w
  # number of nodes
  n <- nrow(S)
  # number of edges
  m <- .5 * n * (n - 1)
  # l1-norm penalty factor
  J <- matrix(1, n, n) / n
  H <- 2 * diag(n) - n * J
  K <- S + alpha * H
  E <- get_incidence_from_adjacency(A)
  R <- t(E) %*% K %*% E
  r <- nrow(R)
  G <- cbind(E, rep(1, n))
  if (record_objective)
    fun <- obj_func(E, K, wk, J)
  # MM-loop
  for (k in c(1:maxiter)) {
    w_aug <- c(wk, 1 / n)
    G_aug_t <- t(G) * w_aug
    G_aug <- t(G_aug_t)
    Q <- G_aug_t %*% solve(G_aug %*% t(G), G_aug)
    Q <- Q[1:m, 1:m]
    wk <- sqrt(diag(Q) / diag(R))
    if (record_objective)
      fun <- c(fun, obj_func(E, K, wk, J))
    werr <- abs(w - wk)
    has_converged <- all(werr <= .5 * reltol * (w + wk)) || all(werr <= abstol)
    if (has_converged && k > 1) break
    w <- wk
  }
  results <- list(Laplacian = E %*% diag(wk) %*% t(E), maxiter = k, converged = has_converged)
  if (record_objective)
    results$objective_function <- fun
  return(results)
}

obj_func <- function(E, K, w, J) {
  n <- ncol(J)
  EWEt <- E %*% diag(w) %*% t(E)
  Gamma <- EWEt + J
  lambda <- eigval_sym(Gamma)[2:n]
  return(sum(E %*% diag(w) %*% t(E)) - sum(log(lambda)))
}
