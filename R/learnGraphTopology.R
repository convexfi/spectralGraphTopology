#' Learn Graph Topology
#'
#' Learns the topology of a K-connected graph given an observed data matrix
#'
#' @param Y T-by-N observed data matrix, where T is the size of each sample
#' and N is the number of samples
#' @param K the number of components of the graph
#' @param w0 initial estimate for the weight vector the graph
#' @param U0 initial estimate for the matrix whose columns are the eigenvectors
#' of the Laplacian matrix
#' @param lamdba0 initial estimate for the vector whose entries are the
#' eigenvalues of the Laplacian matrix
#' @param lb lower bound for the eigenvalues of the Laplacian matrix
#' @param ub upper bound for the eigenvalues of the Laplacian matrix
#' @param alpha tunning parameter
#' @param beta parameter that controls the strength of the regularization term
#' @param rho how much to increase beta after a complete round of iterations
#' @param maxiter the maximum number of iterations for each beta
#' @param w_tol relative tolerance on w
#' @param U_tol relative tolerance on U
#' @param lambda_tol relative tolerance on lambda
#' @param ftol relative tolerance on the objective function
#' @return Theta the learned Laplacian matrix
#' @return fun_seq the objective function value at every iteration
#' @author Convex group - HKUST
#' @references
#' @examples
#' library(learnGraphTopology)
#'
#' # simulate a Laplacian matrix of a single-component graph
#' w <- sample(1:10, 10)
#' Theta <- L(w)
#'
#' # create fake data
#' T <- 10000
#' N <- ncol(Theta)
#' Y <- MASS::mvrnorm(T, rep(0, N), MASS::ginv(Theta))
#'
#' # learn the Laplacian matrix from the simulated data
#' Theta_est <- learnGraphTopology(Y, 1)
#'
#' # show the relative error between the true Laplacian and the learned one
#' norm(Theta - Theta_est, type="F") / norm(Theta, type="F")
#' @export
learnGraphTopology <- function (Y, K, w0 = NA, U0 = NA, lambda0 = NA, lb = 1e-4,
                                ub = 1e4, alpha = 0., beta = .5, rho = .1,
                                maxiter = 5000, w_tol = 1e-4, U_tol = 1e-4,
                                lambda_tol = 1e-4, ftol = 1e-6) {
  N <- ncol(Y)
  T <- nrow(Y)
  S <- crossprod(Y) / T
  H <- alpha * (2. * diag(N) - matrix(1, N, N))
  Km <- S + H

  # define "appropriate" inital guess
  if (any(is.na(w0)))
    w0 <- runif(as.integer(.5 * N * (N - 1)))
  evd <- eigen(L(w0))
  if (any(is.na(U0)))
    U0 <- evd$vectors[, N:1]
  if (any(is.na(lambda0))) {
    lambda0 <- evd$values
    lambda0 <- lambda0[N:1]
    lambda0[1:K] <- 0
  }

  fun0 <- objFunction(w0, U0, lambda0, Km, beta, N, K)
  fun_seq <- c(fun0)

  for (i in 1:maxiter) {
    for (k in 1:maxiter) {
      w <- w_update(w0, U0, beta, lambda0, N, Km)
      U <- U_update(w, N)
      lambda <- lambda_update(lb, ub, beta, U, w, N, K)

      # check tolerance on parameters
      w_err <- norm(w - w0, type="2") / max(1, norm(w, type="2"))
      U_err <- norm(U - U0, type="F") / max(1, norm(U, type="F"))
      lambda_err <- norm(lambda - lambda0, type="2") /
                        max(1, norm(lambda, type="2"))

      if ((w_err < w_tol) & (U_err < U_tol) & (lambda_err < lambda_tol))
        break

      # check tolerance on objective function
      fun <- objFunction(w, U, lambda, Km, beta, N, K)
      ferr <- abs(fun - fun0) / max(1, abs(fun))
      fun_seq <- c(fun, fun_seq)

      if (ferr < ftol)
        break

      fun0 <- fun
      w0 <- w
      U0 <- U
      lambda0 <- lambda
    }
    beta <- beta * (1 + rho)
  }

  return(list(Theta = L(w), fun = fun_seq))
}

objFunction <- function(w, U, lambda, Km, beta, N, K) {
  Theta <- L(w)
  return(sum(-log(lambda[(K+1):N])) + sum(diag(Km %*% Theta)) +
         .5 * beta * norm(Theta - crossprod(t(U) * sqrt(lambda)), type="F")^2)
}
