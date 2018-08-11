spectralGraphTopology <- function (data, K, lb, ub, alpha, beta = .5, pho = .1,
                                   maxiter = 50, w_tol = 1e-2, U_tol = 1e-2,
                                   Lambda_tol = 1e-2, Lw_tol = 1e-2) {
  # Function to solve the graph learning problem with K-components
  #
  # Args:
  #   data: matrix
  #     a n by k matrix, where n is the size of each sample
  #     and k is the number of samples
  #   K: scalar
  #     number of components of the graph
  #   lb, ub: scalars
  #     lower and upper bounds on the eigenvalues of the Laplacian
  #     matrix
  #
  # Returns:
  #   Lw: matrix
  #     the Laplacian matrix
  #
  n <- nrow(data)
  k <- ncol(data)
  S <- data %*% t(data) / k
  H <- alpha * (2. * diag(n) - matrix(1, n, n))
  Km <- S + H

  # define appropriate inital guess
  w1 <- array(1., as.integer(.5 * n * (n - 1)))
  U1 <- diag(n)[, (K+1):n]
  Lambda1 <- array(1., n-K)

  Lw1 <- LOp(w1, n)
  while (iter < maxiter) {
    while (1) {
      w <- w_update(w1, U1, beta, diag(Lambda1), n, Km)
      U <- U_update(w, n)
      Lambda <- Lambda_update(lb, ub, beta, U, w, n, K)

      w_err <- norm(w - w1, type="2") / max(1., norm(w, type=2))
      U_err <- norm(U - U1, type="F") / max(1., norm(U, type="F"))
      Lambda_err <- norm(Lambda - Lambda1, type="2") /
                        max(1., norm(Lambda, type="2"))

      if ((w_err < w_tol) & (U_err < U_tol) & (Lambda_err < Lambda_tol))
        break

      w1 <- w
      U1 <- U
      Lambda1 <- Lambda
    }
    Lw <- LOp(w, n)
    Lw_err <- norm(Lw - Lw1, type="F") / max(1., norm(Lw, type="F"))
    if (Lw_err < Lw_tol)
      break
    Lw1 <- Lw
    iter <- iter + 1
    beta <- beta * (1 + pho)
  }

  return (Lw)
}
