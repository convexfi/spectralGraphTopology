#' Learn Graph Topology
#'
#' Learns the topology of a K-connected graph given an observed data matrix
#'
#' @param Y T-by-N observed data matrix, where T is the size of each sample
#' and N is the number of samples
#' @param K the number of components of the graph
#' @param w0 initial estimate for the weight vector the graph or a string
#' selecting an appropriate method. Available methods are: "qp": solves a
#' simple quadratic program; "naive": sets w0 to the negative of the
#' off-diagonal elements of the generalized precision matrix; "glasso": uses
#' the glasso solution
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
#' @return Lw the learned Laplacian matrix
#' @return W the learned weighted adjacency matrix
#' @return obj_fun the objective function value at every iteration
#' @return loglike the negative loglikelihood function value at every iteration
#' @return w the optimal value of the weight vector
#' @return lambda the optimal value of the eigenvalues
#' @return U the optimal value of the eigenvectors
#'
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
#' T <- 40
#' N <- ncol(Lw)
#' Y <- MASS::mvrnorm(T, rep(0, N), MASS::ginv(Lw))
#'
#' # learn the Laplacian matrix from the simulated data
#' res <- learnGraphTopology(Y, 1)
#'
#' # relative error between the true Laplacian and the learned one
#' norm(Lw - res$Lw, type="F") / norm(Lw, type="F")
#' @export
learnGraphTopology <- function (Y, K, w0 = "qp", lb = 1e-4, ub = 1e4, alpha = 0.,
                                beta = 10., rho = .1, maxiter = 5000, maxiter_beta = 1,
                                w_tol = 1e-6, lambda_tol = 1e-6, U_tol = 1e-6,
                                ftol = 1e-6) {
  # number of samples per node
  T <- nrow(Y)
  # number of nodes
  N <- ncol(Y)
  # sample covariance matrix
  S <- cov(Y)
  # l1-norm penalty factor
  H <- alpha * (2 * diag(N) - matrix(1, N, N))
  Kmat <- S + H
  # find an appropriate inital guess
  Sinv <- MASS::ginv(S)
  if (w0 == "qp") {
    R <- vecLmat(ncol(Sinv))
    qp <- quadprog::solve.QP(crossprod(R), t(R) %*% vec(Sinv), diag(ncol(R)))
    w0 <- qp$solution
  } else if (w0 == "naive") {
    w0 <- pmax(0, Linv(Sinv))
  }
  evd <- eigen(L(w0))
  lambda0 <- pmax(0, evd$values[(N-K):1])
  U0 <- evd$vectors[, (N-K):1]
  # save objective function value at initial guess
  ll0 <- logLikelihood(L(w0), lambda0, Kmat)
  fun0 <- ll0 + logPrior(beta, L(w0), lambda0, U0)
  fun_seq <- c(fun0)
  ll_seq <- c(ll0)
  w_seq <- list(w0)

  for (i in 1:maxiter_beta) {
    for (k in 1:maxiter) {
      w <- w_update(w0, U0, beta, lambda0, N, Kmat)
      U <- U_update(w, N, K)
      lambda <- lambda_update(lb, ub, beta, U, w, N, K)
      # compute relative error on the parameters
      w_err <- norm(w - w0, type="2") / max(1, norm(w, type="2"))
      lambda_err <- norm(lambda - lambda0, type="2") / max(1, norm(lambda, type="2"))
      U_err <- norm(U - U0, type="F") / max(1, norm(U, type="2"))
      # check tolerance on the parameters
      if ((w_err < w_tol) & (lambda_err < lambda_tol) & (U_err < U_tol))
        break
      # compute relative error on the objective function
      ll <- logLikelihood(L(w), lambda, Kmat)
      fun <- ll + logPrior(beta, L(w), lambda, U)
      ferr <- abs(fun - fun0) / max(1, abs(fun))
      # check tolerance on the objective function
      if (ferr < ftol)
        break
      # save and update estimates
      ll_seq <- c(ll_seq, ll)
      fun_seq <- c(fun_seq, fun)
      w_seq <- rlist::list.append(w_seq, w)
      fun0 <- fun
      w0 <- w
      U0 <- U
      lambda0 <- lambda
    }
    beta <- beta * (1 + rho)
  }
  # compute the Laplacian matrix of the best w
  Lw <- L(w)
  # compute the adjancency matrix
  W <- diag(diag(Lw)) - Lw
  return(list(Lw = Lw, W = W, obj_fun = fun_seq, loglike = ll_seq,
              w_seq = w_seq, w = w, lambda = lambda, U = U))
}
