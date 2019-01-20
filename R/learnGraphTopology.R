#' Learn Graph Topology
#'
#' Learns the topology of a K-connected graph given an observed data matrix
#'
#' @param S N-by-N sample covariance matrix, where N is the number of nodes
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
#' @param Lwtol relative tolerance on the Laplacian matrix
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
#' res <- learnGraphTopology(cov(Y))
#'
#' # relative error between the true Laplacian and the learned one
#' norm(Lw - res$Lw, type="F") / norm(Lw, type="F")
#' @export
learnLaplacianGraphTopology <- function(S, K = 1, w0 = "naive", lb = 1e-4, ub = 1e4, alpha = 0.,
                                        beta = 1., beta_max = beta, nbeta = 1, maxiter = 10000,
                                        Lwtol = 1e-4, ftol = 1e-6) {
  # number of nodes
  N <- nrow(S)
  # l1-norm penalty factor
  H <- alpha * (2 * diag(N) - matrix(1, N, N))
  Kmat <- S + H
  # find an appropriate inital guess
  Sinv <- MASS::ginv(S)
  w0 <- w_init(w0, Sinv)
  # compute quantities on the initial guess
  Lw0 <- L(w0)
  #U0 <- eigenvectors(Lw0)[, (K+1):N]
  #lambda0 <- pmin(ub, pmax(lb, c(eigenvalues(Lw0)[(K+1):N])))
  U0 <- laplacian.U_update(Lw0, N, K)
  lambda0 <- laplacian.lambda_update(lb, ub, beta, U0, Lw0, N, K)
  # save objective function value at initial guess
  ll0 <- laplacian.logLikelihood(Lw0, lambda0, Kmat)
  fun0 <- ll0 + laplacian.logPrior(beta, Lw0, lambda0, U0)
  fun_seq <- c(fun0)
  ll_seq <- c(ll0)
  #w_seq <- list(w0)
  time_seq <- c(0)
  pb = txtProgressBar(min = 0, max = maxiter, initial = 0)
  start_time <- proc.time()[3]
  beta_set <- beta * exp(seq(from = 0, to = log(beta_max/beta), length.out = nbeta))
  for (beta in beta_set) {
    for (k in 1:maxiter) {
      setTxtProgressBar(pb, k)
      w <- laplacian.w_update(w0, Lw0, U0, beta, lambda0, N, Kmat)
      Lw <- L(w)
      U <- laplacian.U_update(Lw, N, K)
      lambda <- laplacian.lambda_update(lb, ub, beta, U, Lw, N, K)
      # compute negloglikelihood and objective function values
      ll <- laplacian.logLikelihood(Lw, lambda, Kmat)
      fun <- ll + laplacian.logPrior(beta, Lw, lambda, U)
      # save estimates
      time_seq <- c(time_seq, proc.time()[3] - start_time)
      ll_seq <- c(ll_seq, ll)
      fun_seq <- c(fun_seq, fun)
      #w_seq <- rlist::list.append(w_seq, w)
      # compute the relative error and check the tolerance on the Laplacian
      # matrix and on the objective function
      Lwerr <- norm(Lw - Lw0, type="F") / max(1, norm(Lw0, type="F"))
      ferr <- abs(fun - fun0) / max(1, abs(fun0))
      if ((Lwerr < Lwtol && ferr < ftol) && k > 1)
        break
      # update estimates
      fun0 <- fun
      w0 <- w
      U0 <- U
      lambda0 <- lambda
      Lw0 <- Lw
    }
  }
  # compute the adjancency matrix
  W <- diag(diag(Lw)) - Lw
  return(list(Lw = Lw, W = W, obj_fun = fun_seq, loglike = ll_seq, #w_seq = w_seq,
              w = w, lambda = lambda, U = U, elapsed_time = time_seq,
              convergence = !(k == maxiter), beta = beta))
}


#' @export
learnAdjacencyGraphTopology <- function(S, z = 0, w0 = "naive", alpha = 0., beta = 1.,
                                        beta_max = beta, nbeta = 1, Lips = NULL,
                                        maxiter = 1e4, Awtol = 1e-4, ftol = 1e-6) {
  n <- ncol(S)
  J <- matrix(1/n, n, n)
  # l1-norm penalty factor
  H <- alpha * (2 * diag(n) - matrix(1, n, n))
  Kmat <- S + H
  # compute initial guess
  Sinv <- MASS::ginv(S)
  w0 <- w_init(w0, Sinv)
  if (is.null(Lips))
    Lips <- 1 / prod(eigenvalues(L(w0) + J))
  # compute quantities on the initial guess
  Aw0 <- A(w0)
  V0 <- adjacency.V_update(Aw0, n, z)
  psi0 <- adjacency.psi_update(V0, Aw0)
  # save objective function value at initial guess
  ll0 <- adjacency.logLikelihood(L(w0), Kmat, J)
  fun0 <- ll0 + adjacency.logPrior(beta, Aw0, psi0, V0)
  fun_t <- Inf
  fun_seq <- c(fun0)
  ll_seq <- c(ll0)
  Lips_seq <- c()
  #w_seq <- list(w0)
  time_seq <- c(0)
  pb = txtProgressBar(min = 0, max = maxiter, initial = 0)
  start_time <- proc.time()[3]
  beta_set <- beta * exp(seq(from = 0, to = log(beta_max/beta), length.out = nbeta))
  for (beta in beta_set) {
    for (k in 1:maxiter) {
      setTxtProgressBar(pb, k)
      # we need to make sure that the Lipschitz constant is large enough
      # in order to avoid divergence
      while(1) {
        w <- adjacency.w_update(w0, Aw0, V0, beta, psi0, Kmat, J, Lips)
        fun_t <- tryCatch({
                     adjacency.objectiveFunction(A(w), L(w), V0, psi0, Kmat, J, beta)
                   }, warning = function(warn) return(Inf), error = function(err) return(Inf)
                 )
        Lips_seq <- c(Lips_seq, Lips)
        if (fun_t > fun0) {
          Lips <- 2 * Lips
        } else {
          Lips <- .5 * Lips
          break
        }
      }
      Lw <- L(w)
      Aw <- A(w)
      V <- adjacency.V_update(Aw, n, z)
      psi <- adjacency.psi_update(V, Aw)
      # compute negloglikelihood and objective function values
      ll <- adjacency.logLikelihood(Lw, Kmat, J)
      fun <- ll + adjacency.logPrior(beta, Aw, psi, V)
      # save estimates
      time_seq <- c(time_seq, proc.time()[3] - start_time)
      ll_seq <- c(ll_seq, ll)
      fun_seq <- c(fun_seq, fun)
      #w_seq <- rlist::list.append(w_seq, w)
      # compute the relative error and check the tolerance on the Laplacian
      # matrix and on the objective function
      Awerr <- norm(Aw - Aw0, type="F") / max(1, norm(Aw0, type="F"))
      ferr <- abs(fun - fun0) / max(1, abs(fun0))
      if ((Awerr < Awtol && ferr < ftol) && k > 1)
        break
      # update estimates
      fun0 <- fun
      w0 <- w
      V0 <- V
      psi0 <- psi
      Aw0 <- Aw
    }
  }
  return(list(Aw = Aw, Lw = Lw, obj_fun = fun_seq,
              loglike = ll_seq, w = w, psi = psi, V = V, elapsed_time = time_seq,
              Lips = Lips, Lips_seq = Lips_seq, convergence = !(k == maxiter), beta = beta))
}

