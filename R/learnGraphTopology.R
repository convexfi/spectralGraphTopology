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
#' library(spectralGraphTopology)
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
#' res <- learn_laplacian_matrix(cov(Y))
#'
#' # relative error between the true Laplacian and the learned one
#' norm(Lw - res$Lw, type="F") / norm(Lw, type="F")
#' @export
learn_laplacian_matrix <- function(S, k = 1, w0 = "naive", lb = 0, ub = 1e4, alpha = 0,
                                   beta = 1, maxiter = 10000, Lwtol = 1e-4, ftol = 1e-4,
                                   eig_tol = 1e-9) {
  # number of nodes
  n <- nrow(S)
  # l1-norm penalty factor
  H <- alpha * (2 * diag(n) - matrix(1, n, n))
  K <- S + H
  # find an appropriate inital guess
  Sinv <- MASS::ginv(S)
  w0 <- w_init(w0, Sinv)
  # compute quantities on the initial guess
  Lw0 <- L(w0)
  U0 <- laplacian.U_update(Lw0, k)
  lambda0 <- laplacian.lambda_update(lb, ub, beta, U0, Lw0, k)
  # save objective function value at initial guess
  ll0 <- laplacian.likelihood(Lw0, lambda0, K)
  fun0 <- ll0 + laplacian.prior(beta, Lw0, lambda0, U0)
  fun_seq <- c(fun0)
  ll_seq <- c(ll0)
  beta_seq <- c(beta)
  #w_seq <- list(w0)
  time_seq <- c(0)
  pb = txtProgressBar(min = 0, max = maxiter, initial = 0)
  start_time <- proc.time()[3]
  for (i in 1:maxiter) {
    setTxtProgressBar(pb, i)
    w <- laplacian.w_update(w0, Lw0, U0, beta, lambda0, K)
    Lw <- L(w)
    U <- laplacian.U_update(Lw, k)
    lambda <- laplacian.lambda_update(lb, ub, beta, U, Lw, k)
    # compute negloglikelihood and objective function values
    ll <- laplacian.likelihood(Lw, lambda, K)
    fun <- ll + laplacian.prior(beta, Lw, lambda, U)
    # save estimates
    time_seq <- c(time_seq, proc.time()[3] - start_time)
    ll_seq <- c(ll_seq, ll)
    fun_seq <- c(fun_seq, fun)
    #w_seq <- rlist::list.append(w_seq, w)
    # compute the relative error and check the tolerance on the Laplacian
    # matrix and on the objective function
    Lwerr <- norm(Lw - Lw0, type="F") / max(1, norm(Lw0, type="F"))
    ferr <- abs(fun - fun0) / max(1, abs(fun0))
    eig_vals <- eigenvalues(Lw)
    n_zero_eigenvalues <- sum(abs(eig_vals) < eig_tol)
    if (k < n_zero_eigenvalues)
      beta <- .99 * beta
    else if (k > n_zero_eigenvalues)
      beta <- 1.01 * beta
    if ((Lwerr < Lwtol) && (i > 1) && k == n_zero_eigenvalues)
      break
    # update estimates
    fun0 <- fun
    w0 <- w
    U0 <- U
    lambda0 <- lambda
    Lw0 <- Lw
    beta_seq <- c(beta_seq, beta)
  }
  # compute the adjancency matrix
  Aw <- A(w)
  return(list(Lw = Lw, Aw = Aw, obj_fun = fun_seq, loglike = ll_seq,
              w = w, lambda = lambda, U = U, elapsed_time = time_seq,
              convergence = !(i == maxiter), beta_seq = beta_seq))
}


#' @export
learn_bipartite_graph <- function(S, z = 0, w0 = "naive", alpha = 0., beta = 1e3,
                                  beta_max = beta, nbeta = 1, Lips = NULL,
                                  maxiter = 1e4, Awtol = 1e-4, ftol = 1e-4) {
  # number of nodes
  n <- ncol(S)
  J <- matrix(1/n, n, n)
  # l1-norm penalty factor
  H <- alpha * (2 * diag(n) - matrix(1, n, n))
  K <- S + H
  # compute initial guess
  w0 <- w_init(w0, MASS::ginv(S))
  if (is.null(Lips))
    Lips <- 1 / min(eigenvalues(L(w0) + J))
  # compute quantities on the initial guess
  Aw0 <- A(w0)
  V0 <- bipartite.V_update(Aw0, z)
  psi0 <- bipartite.psi_update(V0, Aw0)
  # save objective function value at initial guess
  #ll0 <- bipartite.loglikelihood(L(w0), K, J)
  #fun0 <- ll0 + bipartite.prior(beta, Aw0, psi0, V0)
  #fun_seq <- c(fun0)
  #ll_seq <- c(ll0)
  fun_seq <- c()
  ll_seq <- c()
  Lips_seq <- c()
  #w_seq <- list(w0)
  time_seq <- c(0)
  pb = txtProgressBar(min = 0, max = maxiter, initial = 0)
  start_time <- proc.time()[3]
  beta_set <- beta * exp(seq(from = 0, to = log(beta_max/beta), length.out = nbeta))
  for (beta in beta_set) {
    ll0 <- bipartite.likelihood(L(w0), K, J)
    fun0 <- ll0 + bipartite.prior(beta, Aw0, psi0, V0)
    fun_seq <- c(fun_seq, fun0)
    ll_seq <- c(ll_seq, ll0)
    for (i in 1:maxiter) {
      setTxtProgressBar(pb, i)
      # we need to make sure that the Lipschitz constant is large enough
      # in order to avoid divergence
      while(1) {
        # compute the update for w
        w <- bipartite.w_update(w0, Aw0, V0, beta, psi0, K, J, Lips)
        # compute the objective function at the updated value of w
        fun_t <- tryCatch({
                     bipartite.obj_fun(A(w), L(w), V0, psi0, K, J, beta)
                   }, warning = function(warn) return(Inf), error = function(err) return(Inf)
                 )
        # check if the current value of the objective function is
        # larger than the previous one
        Lips_seq <- c(Lips_seq, Lips)
        if (fun_t > (fun0 + ftol))
          # in case it is in fact larger, then increase Lips and recompute w
          Lips <- 2 * Lips
        else {
          # otherwise decrease Lips and get outta here!
          Lips <- .5 * Lips
          break
        }
      }
      Lw <- L(w)
      Aw <- A(w)
      V <- bipartite.V_update(Aw, z)
      psi <- bipartite.psi_update(V, Aw)
      # compute negloglikelihood and objective function values
      ll <- bipartite.likelihood(Lw, K, J)
      fun <- ll + bipartite.prior(beta, Aw, psi, V)
      # save measurements of time and objective functions
      time_seq <- c(time_seq, proc.time()[3] - start_time)
      ll_seq <- c(ll_seq, ll)
      fun_seq <- c(fun_seq, fun)
      #w_seq <- rlist::list.append(w_seq, w)
      # compute the relative error and check the tolerance on the Adjacency
      # matrix and on the objective function
      Awerr <- norm(Aw - Aw0, type="F") / max(1, norm(Aw0, type="F"))
      ferr <- abs(fun - fun0) / max(1, abs(fun0))
      if ((Awerr < Awtol && ferr < ftol) && i > 1)
        break
      # update estimates
      fun0 <- fun
      w0 <- w
      V0 <- V
      psi0 <- psi
      Aw0 <- Aw
    }
  }
  return(list(Lw = Lw, Aw = Aw, obj_fun = fun_seq, loglike = ll_seq, w = w,
              psi = psi, V = V, elapsed_time = time_seq, Lips = Lips,
              Lips_seq = Lips_seq, convergence = !(i == maxiter), beta = beta))
}


#' @export
learn_adjacency_and_laplacian <- function(S, z = 0, k = 1, w0 = "naive", alpha = 0.,
                                          beta1 = 1, beta2 = 1e3, lb = 0, ub = 1e4,
                                          maxiter = 1e4, Lwtol = 1e-4, ftol = 1e-4) {
  # number of nodes
  n <- ncol(S)
  # l1-norm penalty factor
  H <- alpha * (2 * diag(n) - matrix(1, n, n))
  K <- S + H
  # compute initial guess
  w0 <- w_init(w0, MASS::ginv(S))
  # compute quantities on the initial guess
  Aw0 <- A(w0)
  Lw0 <- L(w0)
  V0 <- joint.V_update(Aw0, z)
  psi0 <- joint.psi_update(V0, Aw0)
  U0 <- joint.U_update(Lw0, k)
  lambda0 <- joint.lambda_update(lb, ub, beta1, U0, Lw0, k)
  # save objective function value at initial guess
  ll0 <- joint.likelihood(Lw0, lambda0, K)
  fun0 <- ll0 + joint.prior(beta1, beta2, Lw0, Aw0, U0, V0, lambda0, psi0)
  fun_seq <- c(fun0)
  ll_seq <- c(ll0)
  time_seq <- c(0)
  pb = txtProgressBar(min = 0, max = maxiter, initial = 0)
  start_time <- proc.time()[3]
  for (i in c(1:maxiter)) {
    setTxtProgressBar(pb, i)
    w <- joint.w_update(w0, Lw0, Aw0, U0, V0, lambda0, psi0, beta1, beta2, K)
    Lw <- L(w)
    Aw <- A(w)
    U <- joint.U_update(Lw, k)
    V <- joint.V_update(Aw, z)
    lambda <- joint.lambda_update(lb, ub, beta1, U, Lw, k)
    psi <- joint.psi_update(V, Aw)
    # compute negloglikelihood and objective function values
    ll <- joint.likelihood(Lw, lambda, K)
    fun <- ll + joint.prior(beta1, beta2, Lw, Aw, U, V, lambda, psi)
    # save measurements of time and objective functions
    time_seq <- c(time_seq, proc.time()[3] - start_time)
    ll_seq <- c(ll_seq, ll)
    fun_seq <- c(fun_seq, fun)
    #w_seq <- rlist::list.append(w_seq, w)
    # compute the relative error and check the tolerance on the Adjacency
    # matrix and on the objective function
    Lwerr <- norm(Lw - Lw0, type="F") / max(1, norm(Lw0, type="F"))
    ferr <- abs(fun - fun0) / max(1, abs(fun0))
    if ((Lwerr < Lwtol && ferr < ftol) && i > 1)
      break
    # update estimates
    fun0 <- fun
    w0 <- w
    U0 <- U
    V0 <- V
    lambda0 <- lambda
    psi0 <- psi
    Lw0 <- Lw
    Aw0 <- Aw
  }
  return(list(Lw = Lw, Aw = Aw, obj_fun = fun_seq, loglike = ll_seq, w = w,
              psi = psi, lambda = lambda, V = V, U = U, elapsed_time = time_seq,
              convergence = !(i == maxiter)))
}


#' @export
learn_dregular_graph <- function(S, k = 1, w0 = "qp", alpha = 0.,
                                 beta1 = 1, beta2 = 1e3, lb = 0, ub = 1e4,
                                 maxiter = 1e4, Lwtol = 1e-4, ftol = 1e-4) {
  # number of nodes
  n <- ncol(S)
  # l1-norm penalty factor
  H <- alpha * (2 * diag(n) - matrix(1, n, n))
  K <- S + H
  # compute initial guess
  w0 <- w_init(w0, MASS::ginv(S))
  # compute quantities on the initial guess
  Aw0 <- A(w0)
  Lw0 <- L(w0)
  U0 <- dregular.U_update(Lw0, k)
  lambda0 <- dregular.lambda_update(lb, ub, beta1, U0, Lw0, k)
  d0 <- mean(lambda0)
  # save objective function value at initial guess
  ll0 <- dregular.likelihood(Lw0, lambda0, K)
  fun0 <- ll0 + dregular.prior(beta1, beta2, Lw0, Aw0, U0, lambda0, d0)
  fun_seq <- c(fun0)
  ll_seq <- c(ll0)
  time_seq <- c(0)
  pb = txtProgressBar(min = 0, max = maxiter, initial = 0)
  start_time <- proc.time()[3]
  for (i in c(1:maxiter)) {
    setTxtProgressBar(pb, i)
    w <- dregular.w_update(w0, Lw0, Aw0, U0, beta1, beta2, lambda0, d0, K)
    Lw <- L(w)
    Aw <- A(w)
    U <- dregular.U_update(Lw, k)
    lambda <- dregular.lambda_update(lb, ub, beta1 + beta2, U, Lw, k)
    d <- mean(lambda)
    # compute negloglikelihood and objective function values
    ll <- dregular.likelihood(Lw, lambda, K)
    fun <- ll + dregular.prior(beta1, beta2, Lw, Aw, U, lambda, d)
    # save measurements of time and objective functions
    time_seq <- c(time_seq, proc.time()[3] - start_time)
    ll_seq <- c(ll_seq, ll)
    fun_seq <- c(fun_seq, fun)
    #w_seq <- rlist::list.append(w_seq, w)
    # compute the relative error and check the tolerance on the Adjacency
    # matrix and on the objective function
    Lwerr <- norm(Lw - Lw0, type="F") / max(1, norm(Lw0, type="F"))
    ferr <- abs(fun - fun0) / max(1, abs(fun0))
    if ((Lwerr < Lwtol && ferr < ftol) && i > 1)
      break
    # update estimates
    fun0 <- fun
    w0 <- w
    U0 <- U
    lambda0 <- lambda
    Lw0 <- Lw
    Aw0 <- Aw
  }
  return(list(Lw = Lw, Aw = Aw, d = d, obj_fun = fun_seq, loglike = ll_seq, w = w,
              lambda = lambda, U = U, elapsed_time = time_seq,
              convergence = !(i == maxiter)))
}
