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

#' @title Learn the weighted Laplacian matrix of a graph

learn_laplacian_gle_mm <- function(S, A, w0 = "naive", alpha = 0, maxiter = 1000, reltol = 1e-4, abstol = 1e-5) {
  Sinv <- MASS::ginv(S)
  w <- w_init(w0, Sinv)
  wk <- w
  # number of nodes
  n <- nrow(S)
  # number of edges
  m <- .5 * n * (n - 1)
  # l1-norm penalty factor
  H <- alpha * (2 * diag(n) - matrix(1, n, n))
  K <- S + H
  E <- get_incidence_from_adjacency(A)
  print(E)
  R <- t(E) %*% K %*% E
  r <- nrow(R)
  G <- cbind(E, rep(1, n))
  # MM-loop
  for (k in c(1:maxiter)) {
    w_aug <- c(wk, 1 / n)
    G_aug_t <- t(G) * w_aug
    G_aug <- t(G_aug_t)
    Q <- G_aug_t %*% solve(G_aug %*% t(G), G_aug)
    Q <- Q[1:m, 1:m]
    wk <- sqrt(diag(Q) / diag(R))
    werr <- abs(w - wk)
    has_converged <- all(werr <= .5 * reltol * (w + wk)) || all(werr <= abstol)
    if (has_converged) break
    w <- wk
  }
  return(list(Laplacian = E %*% diag(wk) %*% t(E), maxiter = k, converged = has_converged))
}
