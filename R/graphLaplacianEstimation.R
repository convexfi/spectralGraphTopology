# this module contains implementations of the algorithms
# proposed in L. Zhao et. al. "Optimization Algorithms for
# Graph Laplacian Estimation via ADMM and MM", IEEE Trans. Sign. Proc. 2019.

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
  R <- t(E) %*% K %*% E
  r <- nrow(R)
  G <- cbind(E, rep(1, n))
  # MM-loop
  for (k in c(1:maxiter)) {
    w_aug <- c(w, 1 / n)
    G_aug_t <- t(G) * w_aug
    G_aug <- t(G_aug_t)
    Q <- G_aug_t %*% solve(G_aug %*% t(G), G_aug)
    Q <- Q[1:m, 1:m]
    w <- sqrt(diag(Q) * diag(solve(R, rep(1, r))))
    werr <- abs(w - wk)
    has_converged <- all(werr <= .5 * reltol * (w + wk)) || all(werr <= abstol)
    if (has_converged) break
    wk <- w
  }
  return(list(Laplacian = L(w), Adjacency = A(w), w = w))
}
