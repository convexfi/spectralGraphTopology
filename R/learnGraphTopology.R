#' Learn Graph Topology
#'
#' Learns the topology of a K-connected graph given an observed data matrix
#'
#' @param Y T-by-N observed data matrix, where T is the size of each sample
#' and N is the number of samples
#' @param K the number of components of the graph
#' @param w0 initial estimate for the weight vector the graph
#' @param lb lower bound for the eigenvalues of the Laplacian matrix
#' @param ub upper bound for the eigenvalues of the Laplacian matrix
#' @param alpha tunning parameter
#' @param beta parameter that controls the strength of the regularization term
#' @param rho how much to increase beta after a complete round of iterations
#' @param maxiter the maximum number of iterations for each beta
#' @param maxiter_beta the maximum number of iterations for the outer loop
#' @param w_tol relative tolerance on w
#' @param U_tol relative tolerance on U
#' @param lambda_tol relative tolerance on lambda
#' @param ftol relative tolerance on the objective function
#' @return Theta the learned Laplacian matrix
#' @return fun_seq the objective function value at every iteration
#' @return loglike the negative loglikelihood function value at every iteration
#' @return w the optimal value of the weight vector
#' @return lambda the optimal value of the eigenvalues
#' @return U the optimal value of the eigenvectors
#' @author Convex group - HKUST
#' @references our paper soon to be submitted
#' @examples
#' library(learnGraphTopology)
#'
#' # simulate a Laplacian matrix of a single-component graph
#' w <- sample(1:10, 10)
#' Lw <- L(w)
#'
#' # create fake data
#' T <- 10000
#' N <- ncol(Theta)
#' Y <- MASS::mvrnorm(T, rep(0, N), MASS::ginv(Theta))
#'
#' # learn the Laplacian matrix from the simulated data
#' Lw_est <- learnGraphTopology(Y, 1)
#'
#' # show the relative error between the true Laplacian and the learned one
#' norm(Lw - Lw_est, type="F") / norm(Lw, type="F")
#' @export
learnGraphTopology <- function (Y, K, w0 = NA, lb = 1e-4, ub = 1e4, alpha = 0.,
                                beta = .5, rho = .1, maxiter = 5000, maxiter_beta = 1,
                                w_tol = 1e-6, lambda_tol = 1e-6, U_tol = 1e-6,
                                ftol = 1e-6) {
  N <- ncol(Y)
  T <- nrow(Y)
  S <- cov(Y)
  H <- alpha * (2 * diag(N) - matrix(1, N, N))
  Kmat <- S + H

  # find an appropriate inital guess via QP
  if (any(is.na(w0))) {
    Sinv <- MASS::ginv(Kmat)
    R <- vecLmat(ncol(Sinv))
    qp <- quadprog::solve.QP(t(R) %*% R, t(R) %*% vec(Sinv), diag(ncol(R)))
    w0 <- qp$solution
  }
  evd <- eigen(L(w0))
  lambda0 <- pmax(0, evd$values[(N-K):1])
  U0 <- evd$vectors[, (N-K):1]

  ll0 <- logLikelihood(L(w0), lambda0, Kmat)
  fun0 <- ll0 + logPrior(beta, L(w0), lambda0, U0)
  fun_seq <- c(fun0)
  ll_seq <- c(ll0)

  for (i in 1:maxiter_beta) {
    for (k in 1:maxiter) {
      w <- w_update(w0, U0, beta, lambda0, N, Kmat)
      U <- U_update(w, N, K)
      lambda <- lambda_update(lb, ub, beta, U, w, N, K)

      # check tolerance on parameters
      w_err <- norm(w - w0, type="2") / max(1, norm(w, type="2"))
      lambda_err <- norm(lambda - lambda0, type="2") / max(1, norm(lambda, type="2"))
      U_err <- norm(U - U0, type="F") / max(1, norm(U, type="2"))

      if ((w_err < w_tol) & (lambda_err < lambda_tol) & (U_err < U_tol))
        break

      # check tolerance on objective function
      ll <- logLikelihood(L(w), lambda, Kmat)
      fun <- ll + logPrior(beta, L(w), lambda, U)
      ferr <- abs(fun - fun0) / max(1, abs(fun))
      ll_seq <- c(ll_seq, ll)
      fun_seq <- c(fun_seq, fun)

      # check tolerance on objective function
      if (ferr < ftol)
        break

      fun0 <- fun
      w0 <- w
      U0 <- U
      lambda0 <- lambda
    }
    beta <- beta * (1 + rho)
  }
  Lw <- L(w)
  W <- diag(diag(Lw)) - Lw
  return(list(Lw = L(w), W = W, obj_fun = fun_seq, loglike = ll_seq,
              w = w, lambda = lambda, U = U))
}
