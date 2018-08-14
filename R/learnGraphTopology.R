#' Sparse Index Tracking
#'
#' Computes the weights of assets (relative capital allocation) for a sparse approximation of a financial index.
#'
#' @param X m-by-n matrix of net returns (m samples, n assets).
#' @param r m dimensional vector of the net returns of the index.
#' @param lambda sparsity weight factor. Any nonnegative number (suggested range \code{[10^{-8},10^{-6}]}).
#' @param u upper bound of the weights.Default value \code{u <- 1}, i.e., no effective upper bound.
#' @param measure performance measure. Possible values \code{'ete'} (empirical tracking error - default), \code{'dr'} (downside risk),
#' \code{'hete'} (Huber empirical tracking error), and \code{'hdr'} (Huber downside risk).
#' @param hub Huber parameter. Required if \code{measure = 'hete'} or \code{measure = 'hdr'}.
#' @param w0 initial point. If \code{NULL} a uniform allocation is used, i.e., \code{w0 <- rep(1/N, N)}.
#' @param thres threshold value. All the weights less or equal to \code{thres} are set to 0. The default value is \code{1e-9}.
#' @return An n-dimensional vector with allocation weights on the assets.
#' @author Konstantinos Benidis and Daniel P. Palomar
#' @references
#' K. Benidis, Y. Feng, D. P. Palomar, "Sparse Portfolios for High-Dimensional Financial Index Tracking,"
#' \emph{IEEE Transactions on Signal Processing}, vol. 66, no. 1, pp. 155-170, Jan. 2018.
#' @examples
#' library(sparseIndexTracking)
#' library(xts)
#'
#' # load data
#' data(INDEX_2010)
#'
#' # fit portfolio under error measure ETE (Empirical Tracking Error)
#' w_ete <- spIndexTrack(INDEX_2010$X, INDEX_2010$SP500, lambda = 1e-7, u = 0.5, measure = 'ete')
#'
#' # show cardinality achieved
#' cat("Number of assets used:", sum(w_ete > 1e-6))
#'
#' @export
learnGraphTopology <- function (data, K, w1 = NA, U1 = NA, lambda1 = NA,
                                lb = 1e-4, ub = 1e4, alpha = 0.,
                                beta = .5, rho = .1, maxiter = 50,
                                w_tol = 1e-4, U_tol = 1e-4,
                                lambda_tol = 1e-4, Lw_tol = 1e-4) {
  # Solves the graph learning problem with K-components
  #
  # Args:
  #   data: matrix
  #     a N by T matrix, where N is the size of each sample
  #     and T is the number of samples
  #   K: scalar
  #     number of components of the graph
  #   w1: vector
  #     initial estimate for w
  #   U1: matrix
  #     initial estimate for U
  #   lambda1: vector
  #     initial estimate for lambda
  #   lb, ub: scalars
  #     lower and upper bounds on the eigenvalues of the Laplacian
  #     matrix
  #   alpha: scalar
  #   beta: scalar
  #     regularization factor
  #   rho: float
  #     how much to increase beta per iteration, i.e.,
  #     beta <- beta * (1 + rho)
  #   maxiter: integer
  #     maximum number of iterations
  #   w_tol: float
  #     tolerance for the convergence of w
  #   U_tol: float
  #     tolerance for the convergence of U
  #   lambda_tol: float
  #     tolerance for the convergence of lambda
  #   Lw_tol: float
  #     tolerance for the convergence of Lw
  #
  # Returns:
  #   Lw: matrix
  #     the Laplacian matrix
  #
  N <- nrow(data)
  T <- ncol(data)
  S <- tcrossprod(data) / T
  H <- alpha * (2. * diag(N) - matrix(1, N, N))
  Km <- S + H

  # define "appropriate" inital guess
  if (any(is.na(w1)))
    w1 <- array(1., as.integer(.5 * N * (N - 1)))
  Lw1 <- L(w1)
  evd <- eigen(Lw1)
  if (any(is.na(U1)))
    U1 <- evd$vectors[, 1:(N-K)]
  if (any(is.na(lambda1)))
    lambda1 <- evd$values[1:(N-K)]

  iter_1 <- 1
  iter_2 <- 1
  while (iter_1 < maxiter) {
    while (iter_2 < maxiter) {
      w <- w_update(w1, U1, beta, lambda1, N, Km)
      U <- U_update(w, N, K)
      lambda <- lambda_update(lb, ub, beta, U, w, N, K)

      w_err <- norm(w - w1, type="2") / max(1., norm(w, type="2"))
      U_err <- norm(U - U1, type="F") / max(1., norm(U, type="F"))
      lambda_err <- norm(lambda - lambda1, type="2") /
                        max(1., norm(lambda, type="2"))

      if ((w_err < w_tol) & (U_err < U_tol) & (lambda_err < lambda_tol))
        break

      w1 <- w
      U1 <- U
      lambda1 <- lambda
      iter_2 <- iter_2 + 1
    }
    Lw <- L(w)
    Lw_err <- norm(Lw - Lw1, type="F") / max(1., norm(Lw, type="F"))
    if (Lw_err < Lw_tol)
      break
    Lw1 <- Lw
    iter_1 <- iter_1 + 1
    beta <- beta * (1 + rho)
  }

  return (Lw)
}
