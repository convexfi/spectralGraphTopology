#' Learn Graph Topology
#'
#' Learns the topology of a K-connected graph given an observed data matrix
#'
#' @param S either a pxp sample covariance/correlation matrix, or a pxn data
#'        matrix, where p is the number of nodes and n is the number of
#'        features (or data points per node)
#' @param is_data_matrix whether the matrix S should be treated as data matrix
#'        or sample covariance matrix
#' @param m in case is_data_matrix = TRUE, then we build an affinity matrix based
#'        on Nie et. al. 2017, where m is the maximum number of possible connections
#'        for a given node
#' @param k the number of components of the graph
#' @param w0 initial estimate for the weight vector the graph or a string
#'        selecting an appropriate method. Available methods are: "qp": finds w0 that minimizes
#'        ||ginv(S) - L(w0)||_F, w0 >= 0; "naive": takes w0 as the negative of the
#'        off-diagonal elements of the pseudo inverse, setting to 0 any elements s.t.
#'        w0 < 0
#' @param lb lower bound for the eigenvalues of the Laplacian matrix
#' @param ub upper bound for the eigenvalues of the Laplacian matrix
#' @param alpha L1 regularization hyperparameter
#' @param beta regularization hyperparameter for the term ||L(w) - U Lambda U'||^2_F
#' @param beta_max maximum allowed value for beta
#' @param fix_beta whether or not to fix the value of beta. In case this parameter
#'        is set to false, then beta will increase (decrease) depending whether the number of
#'        zero eigenvalues is lesser (greater) than k
#' @param rho how much to increase (decrease) beta in case fix_beta = FALSE
#' @param maxiter the maximum number of iterations
#' @param abstol absolute tolerance on the weight vector w
#' @param reltol relative tolerance on the weight vector w
#' @param eigtol value below which eigenvalues are considered to be zero
#' @param record_objective whether to record the objective function values at
#'        each iteration
#' @param record_weights whether to record the edge values at each iteration
#' @param verbose whether to output a progress bar showing the evolution of the
#'        iterations
#' @return Lw the learned Laplacian matrix
#' @return W the learned weighted adjacency matrix
#' @return obj_fun the objective function value at every iteration
#' @return loglike the negative loglikelihood function value at every iteration
#' @return w the optimal value of the weight vector
#' @return lambda the optimal value of the eigenvalues
#' @return U the optimal value of the eigenvectors
#'
#' @author Ze Vinicius and Daniel Palomar
#' @references our paper soon to be submitted
#' @examples
#' library(spectralGraphTopology)
#' library(clusterSim)
#' library(igraph)
#' set.seed(42)
#' # number of nodes per cluster
#' n <- 50
#' # generate datapoints
#' twomoon <- shapes.two.moon(n)
#' # number of components
#' k <- 2
#' # compute sample correlation matrix
#' S <- crossprod(t(twomoon$data))
#' # estimate underlying graph
#' graph <- learn_laplacian_matrix(S, k = k, beta = .5, verbose = FALSE, abstol = 1e-3)
#' # build network
#' net <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected", weighted = TRUE)
#' # colorify nodes and edges
#' colors <- c("#706FD3", "#FF5252")
#' V(net)$cluster <- twomoon$clusters
#' E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1,
#'                       function(x) ifelse(V(net)$cluster[x[1]] == V(net)$cluster[x[2]],
#'                                         colors[V(net)$cluster[x[1]]], '#000000'))
#' V(net)$color <- colors[twomoon$clusters]
#' # plot nodes
#' plot(net, layout = twomoon$data, vertex.label = NA, vertex.size = 3)
#' @export
learn_laplacian_matrix <- function(S, is_data_matrix = FALSE, k = 1, w0 = "naive", lb = 0, ub = 1e4, alpha = 0,
                                   beta = 1e4, beta_max = 1e6, fix_beta = TRUE, rho = 1e-2, m = 7,
                                   maxiter = 1e4, abstol = 1e-6, reltol = 1e-4, eig_tol = 1e-9,
                                   record_objective = FALSE, record_weights = FALSE, verbose = TRUE) {
  if (is_data_matrix || ncol(S) != nrow(S)) {
    A <- build_initial_graph(S, m = m)
    D <- diag(.5 * colSums(A + t(A)))
    L <- D - .5 * (A + t(A))
    S <- MASS::ginv(L)
    is_data_matrix <- TRUE
  }
  # number of nodes
  n <- nrow(S)
  # l1-norm penalty factor
  H <- alpha * (2 * diag(n) - matrix(1, n, n))
  K <- S + H
  # find an appropriate inital guess
  if (is_data_matrix)
    Sinv <- L
  else
    Sinv <- MASS::ginv(S)
  # if w0 is either "naive" or "qp", compute it, else return w0
  w0 <- w_init(w0, Sinv)
  # compute quantities on the initial guess
  Lw0 <- L(w0)
  U0 <- laplacian.U_update(Lw = Lw0, k = k)
  lambda0 <- laplacian.lambda_update(lb = lb, ub = ub, beta = beta, U = U0,
                                     Lw = Lw0, k = k)
  # save objective function value at initial guess
  if (record_objective) {
    ll0 <- laplacian.likelihood(Lw = Lw0, lambda = lambda0, K = K)
    fun0 <- ll0 + laplacian.prior(beta = beta, Lw = Lw0, lambda = lambda0, U = U0)
    fun_seq <- c(fun0)
    ll_seq <- c(ll0)
  }
  beta_seq <- c(beta)
  if (record_weights)
    w_seq <- list(Matrix::Matrix(w0, sparse = TRUE))
  time_seq <- c(0)
  if (verbose)
    pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta  beta: :beta  kth_eigval: :kth_eigval relerr: :relerr",
                                     total = maxiter, clear = FALSE, width = 120)
  start_time <- proc.time()[3]
  for (i in 1:maxiter) {
    w <- laplacian.w_update(w = w0, Lw = Lw0, U = U0, beta = beta,
                            lambda = lambda0, K = K)
    Lw <- L(w)
    U <- laplacian.U_update(Lw = Lw, k = k)
    lambda <- laplacian.lambda_update(lb = lb, ub = ub, beta = beta, U = U,
                                      Lw = Lw, k = k)
    # compute negloglikelihood and objective function values
    if (record_objective) {
      ll <- laplacian.likelihood(Lw = Lw, lambda = lambda, K = K)
      fun <- ll + laplacian.prior(beta = beta, Lw = Lw, lambda = lambda, U = U)
      ll_seq <- c(ll_seq, ll)
      fun_seq <- c(fun_seq, fun)
    }
    if (record_weights)
      w_seq <- rlist::list.append(w_seq, Matrix::Matrix(w, sparse = TRUE))
    # check for convergence
    werr <- abs(w0 - w)
    has_w_converged <- all(werr <= .5 * reltol * (w + w0)) || all(werr <= abstol)
    time_seq <- c(time_seq, proc.time()[3] - start_time)
    if (verbose || (!fix_beta)) eigvals <- eigval_sym(Lw)
    if (verbose) {
      pb$tick(token = list(beta = beta, kth_eigval = eigvals[k],
                           relerr = 2 * max(werr / (w + w0), na.rm = 'ignore')))
    }
    if (!fix_beta) {
      n_zero_eigenvalues <- sum(abs(eigvals) < eig_tol)
      if (k <= n_zero_eigenvalues)
        beta <- (1 + rho) * beta
      else if (k > n_zero_eigenvalues)
        beta <- beta / (1 + rho)
      if (beta > beta_max)
        beta <- beta_max
      beta_seq <- c(beta_seq, beta)
    }
    if (has_w_converged)
      break
    # update estimates
    w0 <- w
    U0 <- U
    lambda0 <- lambda
    Lw0 <- Lw
  }
  # compute the adjancency matrix
  Aw <- A(w)
  results <- list(Laplacian = Lw, Adjacency = Aw, w = w, lambda = lambda, U = U,
                  Laplacian_eigvals = eigval_sym(Lw), elapsed_time = time_seq,
                  convergence = !(i == maxiter), beta_seq = beta_seq)
  if (record_objective) {
    results$obj_fun <- fun_seq
    results$loglike <- ll_seq
  }
  if (record_weights)
    results$w_seq <- w_seq
  return(results)
}


#' @export
learn_bipartite_graph <- function(S, is_data_matrix = FALSE, z = 0, w0 = "naive", alpha = 0., nu = 1e4,
                                  Lips = NULL, m = 7, rho = 1, maxiter = 1e4, abstol = 1e-6, reltol = 1e-4,
                                  record_weights = FALSE) {
  if (is_data_matrix || ncol(S) != nrow(S)) {
    A <- build_initial_graph(S, m = m)
    D <- diag(.5 * colSums(A + t(A)))
    L <- D - .5 * (A + t(A))
    S <- MASS::ginv(L)
    is_data_matrix <- TRUE
  }
  # number of nodes
  n <- nrow(S)
  # note now that S is always some sort of similarity matrix
  J <- matrix(1/n, n, n)
  # l1-norm penalty factor
  H <- alpha * (2 * diag(n) - matrix(1, n, n))
  K <- S + H
  # compute initial guess
  if (is_data_matrix)
    Sinv <- L
  else
    Sinv <- MASS::ginv(S)
  # if w0 is either "naive" or "qp", compute it, else return w0
  w0 <- w_init(w0, Sinv)
  if (is.null(Lips))
    Lips <- 1 / min(eigval_sym(L(w0) + J))
  # compute quantities on the initial guess
  Aw0 <- A(w0)
  V0 <- bipartite.V_update(Aw0, z)
  psi0 <- bipartite.psi_update(V0, Aw0)
  Lips_seq <- c(Lips)
  time_seq <- c(0)
  start_time <- proc.time()[3]
  ll0 <- bipartite.likelihood(Lw = L(w0), K = K, J = J)
  fun0 <- ll0 + bipartite.prior(nu = nu, Aw = Aw0, psi = psi0, V = V0)
  fun_seq <- c(fun0)
  ll_seq <- c(ll0)
  if (record_weights)
    w_seq <- list(Matrix::Matrix(w0, sparse = TRUE))
  pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta  Lipschitz: :Lips  relerr: :relerr",
                                   total = maxiter, clear = FALSE, width = 100)
  for (i in 1:maxiter) {
    # we need to make sure that the Lipschitz constant is large enough
    # in order to avoid divergence
    while(1) {
      # compute the update for w
      w <- bipartite.w_update(w = w0, Aw = Aw0, V = V0, nu = nu, psi = psi0,
                              K = K, J = J, Lips = Lips)
      # compute the objective function at the updated value of w
      fun_t <- tryCatch({
                   bipartite.obj_fun(Aw = A(w), Lw = L(w), V = V0, psi = psi0,
                                     K = K, J = J, nu = nu)
                 }, warning = function(warn) return(Inf), error = function(err) return(Inf)
               )
      # check if the previous value of the objective function is
      # smaller than the current one
      Lips_seq <- c(Lips_seq, Lips)
      if (fun0 < fun_t)
        # in case it is in fact larger, then increase Lips and recompute w
        Lips <- (1 + rho) * Lips
      else {
        # otherwise decrease Lips and get outta here!
        Lips <- Lips / (1 + rho)
        if (Lips < 1e-12)
          Lips <- 1e-12
        break
      }
    }
    Lw <- L(w)
    Aw <- A(w)
    V <- bipartite.V_update(Aw = Aw, z = z)
    psi <- bipartite.psi_update(V = V, Aw = Aw)
    # compute negloglikelihood and objective function values
    ll <- bipartite.likelihood(Lw = Lw, K = K, J = J)
    fun <- ll + bipartite.prior(nu = nu, Aw = Aw, psi = psi, V = V)
    # save measurements of time and objective functions
    time_seq <- c(time_seq, proc.time()[3] - start_time)
    ll_seq <- c(ll_seq, ll)
    fun_seq <- c(fun_seq, fun)
    # compute the relative error and check the tolerance on the Adjacency
    # matrix and on the objective function
    if (record_weights)
      w_seq <- rlist::list.append(w_seq, Matrix::Matrix(w, sparse = TRUE))
    # check for convergence
    werr <- abs(w0 - w)
    has_w_converged <- (all(werr <= .5 * reltol * (w + w0)) || all(werr <= abstol))
    pb$tick(token = list(Lips = Lips, relerr = 2*max(werr/(w + w0), na.rm = 'ignore')))
    if (has_w_converged)
      break
    # update estimates
    fun0 <- fun
    w0 <- w
    V0 <- V
    psi0 <- psi
    Aw0 <- Aw
  }
  results <- list(Laplacian = Lw, Adjacency = Aw, obj_fun = fun_seq, loglike = ll_seq, w = w,
                  psi = psi, V = V, elapsed_time = time_seq, Lips = Lips,
                  Lips_seq = Lips_seq, convergence = (i < maxiter), nu = nu)
  if (record_weights)
    results$w_seq <- w_seq
  return(results)
}


#' @export
learn_adjacency_and_laplacian <- function(S, is_data_matrix = FALSE, z = 0, k = 1,
                                          w0 = "naive", m = 7, alpha = 0., beta = 1e4,
                                          rho = 1e-2, fix_beta = TRUE, beta_max = 1e6, nu = 1e4,
                                          lb = 0, ub = 1e4, maxiter = 1e4, abstol = 1e-6,
                                          reltol = 1e-4, eig_tol = 1e-9,
                                          record_weights = FALSE, record_objective = FALSE) {
  if (is_data_matrix || ncol(S) != nrow(S)) {
    A <- build_initial_graph(S, m = m)
    D <- diag(.5 * colSums(A + t(A)))
    L <- D - .5 * (A + t(A))
    S <- MASS::ginv(L)
    is_data_matrix <- TRUE
  }
  # number of nodes
  n <- nrow(S)
  # l1-norm penalty factor
  H <- alpha * (2 * diag(n) - matrix(1, n, n))
  K <- S + H
  # find an appropriate inital guess
  if (is_data_matrix)
    Sinv <- L
  else
    Sinv <- MASS::ginv(S)
  w0 <- w_init(w0, Sinv)
  # compute quantities on the initial guess
  Aw0 <- A(w0)
  Lw0 <- L(w0)
  V0 <- joint.V_update(Aw0, z)
  psi0 <- joint.psi_update(V0, Aw0)
  U0 <- joint.U_update(Lw0, k)
  lambda0 <- joint.lambda_update(lb, ub, beta, U0, Lw0, k)
  if (record_objective) {
    # save objective function value at initial guess
    ll0 <- joint.likelihood(Lw0, lambda0, K)
    fun0 <- ll0 + joint.prior(beta, nu, Lw0, Aw0, U0, V0, lambda0, psi0)
    fun_seq <- c(fun0)
    ll_seq <- c(ll0)
  }
  beta_seq <- c(beta)
  time_seq <- c(0)
  start_time <- proc.time()[3]
  if (record_weights)
    w_seq <- list(Matrix::Matrix(w0, sparse = TRUE))
  pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta  beta: :beta  kth_eigval: :kth_eigval  relerr: :relerr",
                                   total = maxiter, clear = FALSE, width = 120)
  for (i in c(1:maxiter)) {
    w <- joint.w_update(w0, Lw0, Aw0, U0, V0, lambda0, psi0, beta, nu, K)
    Lw <- L(w)
    Aw <- A(w)
    U <- joint.U_update(Lw, k)
    V <- joint.V_update(Aw, z)
    lambda <- joint.lambda_update(lb, ub, beta, U, Lw, k)
    psi <- joint.psi_update(V, Aw)
    time_seq <- c(time_seq, proc.time()[3] - start_time)
    if (record_objective) {
      ll <- joint.likelihood(Lw, lambda, K)
      fun <- ll + joint.prior(beta, nu, Lw, Aw, U, V, lambda, psi)
      ll_seq <- c(ll_seq, ll)
      fun_seq <- c(fun_seq, fun)
    }
    if (record_weights)
      w_seq <- rlist::list.append(w_seq, Matrix::Matrix(w, sparse = TRUE))
    werr <- abs(w0 - w)
    has_w_converged <- (all(werr <= .5 * reltol * (w + w0)) || all(werr <= abstol))
    time_seq <- c(time_seq, proc.time()[3] - start_time)
    eigvals <- eigval_sym(Lw)
    pb$tick(token = list(beta = beta, kth_eigval = eigvals[k], relerr = 2*max(werr/(w + w0), na.rm = 'ignore')))
    if (!fix_beta) {
      n_zero_eigenvalues <- sum(abs(eigvals) < eig_tol)
      if (k < n_zero_eigenvalues)
        beta <- (1 + rho) * beta
      else if (k > n_zero_eigenvalues)
        beta <- beta / (1 + rho)
      if (beta > beta_max)
        beta <- beta_max
      beta_seq <- c(beta_seq, beta)
    }
    if (has_w_converged)
      break
    # update estimates
    w0 <- w
    U0 <- U
    V0 <- V
    lambda0 <- lambda
    psi0 <- psi
    Lw0 <- Lw
    Aw0 <- Aw
  }
  results <- list(Laplacian = Lw, Adjacency = Aw, w = w, psi = psi,
                  lambda = lambda, V = V, U = U, elapsed_time = time_seq,
                  beta_seq = beta_seq, nu = nu, convergence = (i < maxiter))
  if (record_objective) {
    results$obj_fun <- fun_seq
    results$loglike <- ll_seq
  }
  if (record_weights)
    results$w_seq <- w_seq
  return(results)
}


#' @export
learn_dregular_graph <- function(S, is_data_matrix = FALSE, d0 = NULL, k = 1, w0 = "qp",
                                 alpha = 0, beta = 1e3, beta_max = 1e6, rho = 1e-2,
                                 fix_beta = FALSE, eta = 1e3, lb = 0, ub = 1e4,
                                 maxiter = 1e4, abstol = 1e-6, reltol = 1e-4,
                                 eig_tol = 1e-9, record_objective = FALSE,
                                 record_weights = FALSE) {
  if (is_data_matrix || ncol(S) != nrow(S)) {
    A <- build_initial_graph(S, m = m)
    D <- diag(.5 * colSums(A + t(A)))
    L <- D - .5 * (A + t(A))
    S <- MASS::ginv(L)
    is_data_matrix <- TRUE
  }
  # number of nodes
  n <- nrow(S)
  # l1-norm penalty factor
  H <- alpha * (2 * diag(n) - matrix(1, n, n))
  K <- S + H
  # compute initial guess
  if (is_data_matrix)
    Sinv <- L
  else
    Sinv <- MASS::ginv(S)
  # if w0 is either "naive" or "qp", compute it, else return w0
  w0 <- w_init(w0, Sinv)
  # compute quantities on the initial guess
  Aw0 <- A(w0)
  Lw0 <- L(w0)
  U0 <- dregular.U_update(Lw = Lw0, k = k)
  lambda0 <- dregular.lambda_update(lb = lb, ub = ub, beta = beta, U = U0,
                                    Lw = Lw0, k = k)
  if (is.null(d0))
    d0 <- mean(diag(Lw0))
  # save objective function value at initial guess
  if (record_objective) {
    ll0 <- dregular.likelihood(Lw = Lw0, lambda = lambda0, K = K)
    fun0 <- ll0 + dregular.prior(beta = beta, eta = eta, Lw = Lw0, Aw = Aw0,
                                 U = U0, lambda = lambda0, d = d0)
    fun_seq <- c(fun0)
    ll_seq <- c(ll0)
  }
  time_seq <- c(0)
  start_time <- proc.time()[3]
  beta_seq <- c(beta)
  if (record_weights)
    w_seq <- list(Matrix::Matrix(w0, sparse = TRUE))
  pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta  beta: :beta  kth_eigval: :kth_eigval  relerr: :relerr",
                                   total = maxiter, clear = FALSE, width = 100)
  for (i in c(1:maxiter)) {
    w <- dregular.w_update(w = w0, Lw = Lw0, Aw = Aw0, U = U0, beta = beta,
                           eta = eta, lambda = lambda0, d = d0, K = K)
    Lw <- L(w)
    Aw <- A(w)
    U <- dregular.U_update(Lw = Lw, k = k)
    lambda <- dregular.lambda_update(lb = lb, ub = ub, beta = beta, U = U,
                                     Lw = Lw, k = k)
    d <- mean(diag(Lw))
    if (record_objective) {
      # compute negloglikelihood and objective function values
      ll <- dregular.likelihood(Lw, lambda, K)
      fun <- ll + dregular.prior(beta, eta, Lw, Aw, U, lambda, d)
      # save measurements of time and objective functions
      ll_seq <- c(ll_seq, ll)
      fun_seq <- c(fun_seq, fun)
    }
    if (record_weights)
      w_seq <- rlist::list.append(w_seq, Matrix::Matrix(w, sparse = TRUE))
    time_seq <- c(time_seq, proc.time()[3] - start_time)
    eigvals <- eigval_sym(Lw)
    if (!fix_beta) {
      n_zero_eigenvalues <- sum(abs(eigvals) < eig_tol)
      if (k < n_zero_eigenvalues)
        beta <- (1 + rho) * beta
      else if (k > n_zero_eigenvalues)
        beta <- beta / (1 + rho)
      if (beta > beta_max)
        beta <- beta_max
      beta_seq <- c(beta_seq, beta)
    }
    werr <- abs(w0 - w)
    has_w_converged <- (all(werr <= .5 * reltol * (w + w0)) || all(werr <= abstol))
    pb$tick(token = list(beta = beta, kth_eigval = eigvals[k],
                         relerr = 2 * max(werr / (w + w0), na.rm = 'ignore')))
    if (has_w_converged)
      break
    # update estimates
    d0 <- d
    w0 <- w
    U0 <- U
    lambda0 <- lambda
    Lw0 <- Lw
    Aw0 <- Aw
  }
  results <- list(Laplacian = Lw, Adjacency = Aw, w = w, lambda = lambda,
                  U = U, d = d, elapsed_time = time_seq, beta_seq = beta_seq,
                  convergence = (i < maxiter))
  if (record_objective) {
    results$obj_fun <- fun_seq
    results$loglike <- ll_seq
  }
  if (record_weights)
    results$w_seq <- w_seq
  return(results)
}
