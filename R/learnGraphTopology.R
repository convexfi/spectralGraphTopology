#' @title Learn the Laplacian matrix of a k-component graph
#'
#' Learns a k-component graph on the basis of an observed data matrix.
#' Check out https://mirca.github.io/spectralGraphTopology for code examples.
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
#' @param alpha reweighted l1-norm regularization hyperparameter
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
#' @return A list containing possibly the following elements:
#' \item{\code{laplacian}}{the estimated Laplacian Matrix}
#' \item{\code{adjacency}}{the estimated Adjacency Matrix}
#' \item{\code{w}}{the estimated weight vector}
#' \item{\code{lambda}}{optimization variable accounting for the eigenvalues of the Laplacian matrix}
#' \item{\code{U}}{eigenvectors of the estimated Laplacian matrix}
#' \item{\code{elapsed_time}}{elapsed time recorded at every iteration}
#' \item{\code{beta_seq}}{sequence of values taken by beta in case fix_beta = FALSE}
#' \item{\code{convergence}}{boolean flag to indicate whether or not the optimization converged}
#' \item{\code{obj_fun}}{values of the objective function at every iteration in case record_objective = TRUE}
#' \item{\code{negloglike}}{values of the negative loglikelihood at every iteration in case record_objective = TRUE}
#' \item{\code{w_seq}}{sequence of weight vectors at every iteration in case record_weights = TRUE}
#' @author Ze Vinicius and Daniel Palomar
#' @references S. Kumar, J. Ying, J. V. M. Cardoso, D. P. Palomar. A unified
#'             framework for structured graph learning via spectral constraints.
#'             Journal of Machine Learning Research, 2020.
#'             http://jmlr.org/papers/v21/19-276.html
#' @examples
#' # design true Laplacian
#' Laplacian <- rbind(c(1, -1, 0, 0),
#'                    c(-1, 1, 0, 0),
#'                    c(0, 0, 1, -1),
#'                    c(0, 0, -1, 1))
#' n <- ncol(Laplacian)
#' # sample data from multivariate Gaussian
#' Y <- MASS::mvrnorm(n * 500, rep(0, n), MASS::ginv(Laplacian))
#' # estimate graph on the basis of sampled data
#' graph <- learn_k_component_graph(cov(Y), k = 2, beta = 10)
#' graph$laplacian

#' @export
learn_k_component_graph <- function(S, is_data_matrix = FALSE, k = 1, w0 = "naive", lb = 0, ub = 1e4, alpha = 0,
                                    beta = 1e4, beta_max = 1e6, fix_beta = TRUE, rho = 1e-2, m = 7, eps = 1e-4,
                                    maxiter = 1e4, abstol = 1e-6, reltol = 1e-4, eigtol = 1e-9,
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
  # find an appropriate inital guess
  if (is_data_matrix)
    Sinv <- L
  else
    Sinv <- MASS::ginv(S)
  # if w0 is either "naive" or "qp", compute it, else return w0
  w0 <- w_init(w0, Sinv)
  # compute quantities on the initial guess
  Lw0 <- L(w0)
  # l1-norm penalty factor
  # H <- alpha * (2 * diag(n) - matrix(1, n, n))
  H <- alpha * (diag(n) - matrix(1, n, n))
  K <- S + H
  U0 <- laplacian.U_update(Lw = Lw0, k = k)
  lambda0 <- laplacian.lambda_update(lb = lb, ub = ub, beta = beta, U = U0,
                                     Lw = Lw0, k = k)
  # save objective function value at initial guess
  if (record_objective) {
    ll0 <- laplacian.negloglikelihood(Lw = Lw0, lambda = lambda0, K = K)
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
      ll <- laplacian.negloglikelihood(Lw = Lw, lambda = lambda, K = K)
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
      n_zero_eigenvalues <- sum(abs(eigvals) < eigtol)
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
    K <- S + H / (-Lw + eps)
  }
  # compute the adjancency matrix
  Aw <- A(w)
  results <- list(laplacian = Lw, adjacency = Aw, w = w, lambda = lambda, U = U,
                  elapsed_time = time_seq, convergence = has_w_converged,
                  beta_seq = beta_seq)
  if (record_objective) {
    results$obj_fun <- fun_seq
    results$negloglike <- ll_seq
  }
  if (record_weights)
    results$w_seq <- w_seq
  return(results)
}


learn_cospectral_graph <- function(S, lambda, k = 1, is_data_matrix = FALSE, w0 = "naive", alpha = 0,
                                   beta = 1e4, beta_max = 1e6, fix_beta = TRUE, rho = 1e-2, m = 7,
                                   maxiter = 1e4, abstol = 1e-6, reltol = 1e-4, eigtol = 1e-9,
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
  # save objective function value at initial guess
  if (record_objective) {
    ll0 <- laplacian.negloglikelihood(Lw = Lw0, lambda = lambda, K = K)
    fun0 <- ll0 + laplacian.prior(beta = beta, Lw = Lw0, lambda = lambda, U = U0)
    fun_seq <- c(fun0)
    ll_seq <- c(ll0)
  }
  beta_seq <- c(beta)
  if (record_weights)
    w_seq <- list(Matrix::Matrix(w0, sparse = TRUE))
  time_seq <- c(0)
  if (verbose)
    pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta  beta: :beta  relerr: :relerr",
                                     total = maxiter, clear = FALSE, width = 120)
  start_time <- proc.time()[3]
  for (i in 1:maxiter) {
    w <- laplacian.w_update(w = w0, Lw = Lw0, U = U0, beta = beta,
                            lambda = lambda, K = K)
    Lw <- L(w)
    U <- laplacian.U_update(Lw = Lw, k = k)
    # compute negloglikelihood and objective function values
    if (record_objective) {
      ll <- laplacian.negloglikelihood(Lw = Lw, lambda = lambda, K = K)
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
      n_zero_eigenvalues <- sum(abs(eigvals) < eigtol)
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
    Lw0 <- Lw
  }
  # compute the adjancency matrix
  Aw <- A(w)
  results <- list(laplacian = Lw, adjacency = Aw, w = w, lambda = lambda, U = U,
                  elapsed_time = time_seq, convergence = has_w_converged,
                  beta_seq = beta_seq)
  if (record_objective) {
    results$obj_fun <- fun_seq
    results$negloglike <- ll_seq
  }
  if (record_weights)
    results$w_seq <- w_seq
  return(results)
}


#' @title Learn a bipartite graph
#'
#' Learns a bipartite graph on the basis of an observed data matrix
#'
#' @param S either a pxp sample covariance/correlation matrix, or a pxn data
#'        matrix, where p is the number of nodes and n is the number of
#'        features (or data points per node)
#' @param is_data_matrix whether the matrix S should be treated as data matrix
#'        or sample covariance matrix
#' @param z the number of zero eigenvalues for the Adjancecy matrix
#' @param w0 initial estimate for the weight vector the graph or a string
#'        selecting an appropriate method. Available methods are: "qp": finds w0 that minimizes
#'        ||ginv(S) - L(w0)||_F, w0 >= 0; "naive": takes w0 as the negative of the
#'        off-diagonal elements of the pseudo inverse, setting to 0 any elements s.t.
#'        w0 < 0
#' @param alpha L1 regularization hyperparameter
#' @param m in case is_data_matrix = TRUE, then we build an affinity matrix based
#'        on Nie et. al. 2017, where m is the maximum number of possible connections
#'        for a given node
#' @param nu regularization hyperparameter for the term ||A(w) - V Psi V'||^2_F
#' @param maxiter the maximum number of iterations
#' @param abstol absolute tolerance on the weight vector w
#' @param reltol relative tolerance on the weight vector w
#' @param record_weights whether to record the edge values at each iteration
#' @param verbose whether to output a progress bar showing the evolution of the
#'        iterations
#' @return A list containing possibly the following elements:
#' \item{\code{laplacian}}{the estimated Laplacian Matrix}
#' \item{\code{adjacency}}{the estimated Adjacency Matrix}
#' \item{\code{w}}{the estimated weight vector}
#' \item{\code{psi}}{optimization variable accounting for the eigenvalues of the Adjacency matrix}
#' \item{\code{V}}{eigenvectors of the estimated Adjacency matrix}
#' \item{\code{elapsed_time}}{elapsed time recorded at every iteration}
#' \item{\code{convergence}}{boolean flag to indicate whether or not the optimization converged}
#' \item{\code{obj_fun}}{values of the objective function at every iteration in case record_objective = TRUE}
#' \item{\code{negloglike}}{values of the negative loglikelihood at every iteration in case record_objective = TRUE}
#' \item{\code{w_seq}}{sequence of weight vectors at every iteration in case record_weights = TRUE}
#' @author Ze Vinicius and Daniel Palomar
#' @references S. Kumar, J. Ying, J. V. M. Cardoso, D. P. Palomar. A unified
#'             framework for structured graph learning via spectral constraints.
#'             Journal of Machine Learning Research, 2020.
#'             http://jmlr.org/papers/v21/19-276.html
#' @examples
#' library(spectralGraphTopology)
#' library(igraph)
#' library(viridis)
#' library(corrplot)
#' set.seed(42)
#' n1 <- 10
#' n2 <- 6
#' n <- n1 + n2
#' pc <- .9
#' bipartite <- sample_bipartite(n1, n2, type="Gnp", p = pc, directed=FALSE)
#' # randomly assign edge weights to connected nodes
#' E(bipartite)$weight <- runif(gsize(bipartite), min = 0, max = 1)
#' # get true Laplacian and Adjacency
#' Ltrue <- as.matrix(laplacian_matrix(bipartite))
#' Atrue <- diag(diag(Ltrue)) - Ltrue
#' # get samples
#' Y <- MASS::mvrnorm(100 * n, rep(0, n), Sigma = MASS::ginv(Ltrue))
#' # compute sample covariance matrix
#' S <- cov(Y)
#' # estimate Adjacency matrix
#' graph <- learn_bipartite_graph(S, z = 4, verbose = FALSE)
#' graph$Adjacency[graph$Adjacency < 1e-3] <- 0
#' # Plot Adjacency matrices: true, noisy, and estimated
#' corrplot(Atrue / max(Atrue), is.corr = FALSE, method = "square",
#'          addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
#' corrplot(graph$Adjacency / max(graph$Adjacency), is.corr = FALSE,
#'          method = "square", addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
#' # build networks
#' estimated_bipartite <- graph_from_adjacency_matrix(graph$Adjacency,
#'                                                    mode = "undirected",
#'                                                    weighted = TRUE)
#' V(estimated_bipartite)$type <- c(rep(0, 10), rep(1, 6))
#' la = layout_as_bipartite(estimated_bipartite)
#' colors <- viridis(20, begin = 0, end = 1, direction = -1)
#' c_scale <- colorRamp(colors)
#' E(estimated_bipartite)$color = apply(
#'   c_scale(E(estimated_bipartite)$weight / max(E(estimated_bipartite)$weight)), 1,
#'                           function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
#' E(bipartite)$color = apply(c_scale(E(bipartite)$weight / max(E(bipartite)$weight)), 1,
#'                       function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
#' la = la[, c(2, 1)]
#' # Plot networks: true and estimated
#' plot(bipartite, layout = la, vertex.color=c("red","black")[V(bipartite)$type + 1],
#'      vertex.shape = c("square", "circle")[V(bipartite)$type + 1],
#'      vertex.label = NA, vertex.size = 5)
#' plot(estimated_bipartite, layout = la,
#'      vertex.color=c("red","black")[V(estimated_bipartite)$type + 1],
#'      vertex.shape = c("square", "circle")[V(estimated_bipartite)$type + 1],
#'      vertex.label = NA, vertex.size = 5)
#' @export
learn_bipartite_graph <- function(S, is_data_matrix = FALSE, z = 0, nu = 1e4, alpha = 0.,
                                  w0 = "naive", m = 7, maxiter = 1e4, abstol = 1e-6, reltol = 1e-4,
                                  record_weights = FALSE, verbose = TRUE) {
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
  w0 <- w_init(w0, Sinv) + 1e-4
  Lips <- 1 / min(eigval_sym(L(w0) + J))
  # compute quantities on the initial guess
  Aw0 <- A(w0)
  V0 <- bipartite.V_update(Aw0, z)
  psi0 <- bipartite.psi_update(V0, Aw0)
  Lips_seq <- c(Lips)
  time_seq <- c(0)
  start_time <- proc.time()[3]
  ll0 <- bipartite.negloglikelihood(Lw = L(w0), K = K, J = J)
  fun0 <- ll0 + bipartite.prior(nu = nu, Aw = Aw0, psi = psi0, V = V0)
  fun_seq <- c(fun0)
  ll_seq <- c(ll0)
  if (record_weights)
    w_seq <- list(Matrix::Matrix(w0, sparse = TRUE))
  if (verbose)
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
      Lw <- L(w)
      fun_t <- tryCatch({
                   bipartite.obj_fun(Aw = A(w), Lw = Lw, V = V0, psi = psi0,
                                     K = K, J = J, nu = nu)
                 }, warning = function(warn) return(Inf), error = function(err) return(Inf)
               )
      chol_status <- try(chol(Lw + J), silent = TRUE)
      chol_error <- ifelse(class(chol_status) == "try-error", TRUE, FALSE)
      if (chol_error[1])
        is_disconnected <- TRUE
      else
        is_disconnected <- FALSE
      # check if the previous value of the objective function is
      # smaller than the current one
      Lips_seq <- c(Lips_seq, Lips)
      if (fun0 < fun_t | is_disconnected)
        # in case it is in fact larger, then increase Lips and recompute w
        Lips <- 2 * Lips
      else {
        # otherwise decrease Lips and get outta here!
        Lips <- .5 * Lips
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
    ll <- bipartite.negloglikelihood(Lw = Lw, K = K, J = J)
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
    if (verbose)
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
  results <- list(laplacian = Lw, adjacency = Aw, obj_fun = fun_seq, negloglike = ll_seq, w = w,
                  psi = psi, V = V, elapsed_time = time_seq, Lips = Lips,
                  Lips_seq = Lips_seq, convergence = (i < maxiter), nu = nu)
  if (record_weights)
    results$w_seq <- w_seq
  return(results)
}


#' @title Learns a bipartite k-component graph
#'
#' Jointly learns the Laplacian and Adjacency matrices of a graph on the basis
#' of an observed data matrix
#'
#' @param S either a pxp sample covariance/correlation matrix, or a pxn data
#'        matrix, where p is the number of nodes and n is the number of
#'        features (or data points per node)
#' @param is_data_matrix whether the matrix S should be treated as data matrix
#'        or sample covariance matrix
#' @param z the number of zero eigenvalues for the Adjancecy matrix
#' @param k the number of components of the graph
#' @param w0 initial estimate for the weight vector the graph or a string
#'        selecting an appropriate method. Available methods are: "qp": finds w0 that minimizes
#'        ||ginv(S) - L(w0)||_F, w0 >= 0; "naive": takes w0 as the negative of the
#'        off-diagonal elements of the pseudo inverse, setting to 0 any elements s.t.
#'        w0 < 0
#' @param m in case is_data_matrix = TRUE, then we build an affinity matrix based
#'        on Nie et. al. 2017, where m is the maximum number of possible connections
#'        for a given node
#' @param alpha L1 regularization hyperparameter
#' @param beta regularization hyperparameter for the term ||L(w) - U Lambda U'||^2_F
#' @param rho how much to increase (decrease) beta in case fix_beta = FALSE
#' @param fix_beta whether or not to fix the value of beta. In case this parameter
#'        is set to false, then beta will increase (decrease) depending whether the number of
#'        zero eigenvalues is lesser (greater) than k
#' @param beta_max maximum allowed value for beta
#' @param nu regularization hyperparameter for the term ||A(w) - V Psi V'||^2_F
#' @param lb lower bound for the eigenvalues of the Laplacian matrix
#' @param ub upper bound for the eigenvalues of the Laplacian matrix
#' @param maxiter the maximum number of iterations
#' @param abstol absolute tolerance on the weight vector w
#' @param reltol relative tolerance on the weight vector w
#' @param eigtol value below which eigenvalues are considered to be zero
#' @param record_objective whether to record the objective function values at
#'        each iteration
#' @param record_weights whether to record the edge values at each iteration
#' @param verbose whether to output a progress bar showing the evolution of the
#'        iterations
#'
#' @return A list containing possibly the following elements:
#' \item{\code{laplacian}}{the estimated Laplacian Matrix}
#' \item{\code{adjacency}}{the estimated Adjacency Matrix}
#' \item{\code{w}}{the estimated weight vector}
#' \item{\code{psi}}{optimization variable accounting for the eigenvalues of the Adjacency matrix}
#' \item{\code{lambda}}{optimization variable accounting for the eigenvalues of the Laplacian matrix}
#' \item{\code{V}}{eigenvectors of the estimated Adjacency matrix}
#' \item{\code{U}}{eigenvectors of the estimated Laplacian matrix}
#' \item{\code{elapsed_time}}{elapsed time recorded at every iteration}
#' \item{\code{beta_seq}}{sequence of values taken by beta in case fix_beta = FALSE}
#' \item{\code{convergence}}{boolean flag to indicate whether or not the optimization converged}
#' \item{\code{obj_fun}}{values of the objective function at every iteration in case record_objective = TRUE}
#' \item{\code{negloglike}}{values of the negative loglikelihood at every iteration in case record_objective = TRUE}
#' \item{\code{w_seq}}{sequence of weight vectors at every iteration in case record_weights = TRUE}
#' @author Ze Vinicius and Daniel Palomar
#' @references S. Kumar, J. Ying, J. V. M. Cardoso, D. P. Palomar. A unified
#'             framework for structured graph learning via spectral constraints.
#'             Journal of Machine Learning Research, 2020.
#'             http://jmlr.org/papers/v21/19-276.html
#' @examples
#' library(spectralGraphTopology)
#' library(igraph)
#' library(viridis)
#' library(corrplot)
#' set.seed(42)
#' w <- c(1, 0, 0, 1, 0, 1) * runif(6)
#' Laplacian <- block_diag(L(w), L(w))
#' Atrue <- diag(diag(Laplacian)) - Laplacian
#' bipartite <- graph_from_adjacency_matrix(Atrue, mode = "undirected", weighted = TRUE)
#' n <- ncol(Laplacian)
#' Y <- MASS::mvrnorm(40 * n, rep(0, n), MASS::ginv(Laplacian))
#' graph <- learn_bipartite_k_component_graph(cov(Y), k = 2, beta = 1e2, nu = 1e2, verbose = FALSE)
#' graph$Adjacency[graph$Adjacency < 1e-2] <- 0
#' # Plot Adjacency matrices: true, noisy, and estimated
#' corrplot(Atrue / max(Atrue), is.corr = FALSE, method = "square", addgrid.col = NA, tl.pos = "n",
#'          cl.cex = 1.25)
#' corrplot(graph$Adjacency / max(graph$Adjacency), is.corr = FALSE, method = "square",
#'          addgrid.col = NA, tl.pos = "n", cl.cex = 1.25)
#' # Plot networks
#' estimated_bipartite <- graph_from_adjacency_matrix(graph$Adjacency, mode = "undirected",
#'                                                    weighted = TRUE)
#' V(bipartite)$type <- rep(c(TRUE, FALSE), 4)
#' V(estimated_bipartite)$type <- rep(c(TRUE, FALSE), 4)
#' la = layout_as_bipartite(estimated_bipartite)
#' colors <- viridis(20, begin = 0, end = 1, direction = -1)
#' c_scale <- colorRamp(colors)
#' E(estimated_bipartite)$color = apply(
#'                c_scale(E(estimated_bipartite)$weight / max(E(estimated_bipartite)$weight)), 1,
#'                                      function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
#' E(bipartite)$color = apply(c_scale(E(bipartite)$weight / max(E(bipartite)$weight)), 1,
#'                            function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
#' la = la[, c(2, 1)]
#' # Plot networks: true and estimated
#' plot(bipartite, layout = la,
#'      vertex.color = c("red","black")[V(bipartite)$type + 1],
#'      vertex.shape = c("square", "circle")[V(bipartite)$type + 1],
#'      vertex.label = NA, vertex.size = 5)
#' plot(estimated_bipartite, layout = la,
#'      vertex.color = c("red","black")[V(estimated_bipartite)$type + 1],
#'      vertex.shape = c("square", "circle")[V(estimated_bipartite)$type + 1],
#'      vertex.label = NA, vertex.size = 5)

#' @export
learn_bipartite_k_component_graph <- function(S, is_data_matrix = FALSE, z = 0, k = 1,
                                              w0 = "naive", m = 7, alpha = 0., beta = 1e4,
                                              rho = 1e-2, fix_beta = TRUE, beta_max = 1e6, nu = 1e4,
                                              lb = 0, ub = 1e4, maxiter = 1e4, abstol = 1e-6,
                                              reltol = 1e-4, eigtol = 1e-9,
                                              record_weights = FALSE, record_objective = FALSE, verbose = TRUE) {
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
    ll0 <- joint.negloglikelihood(Lw0, lambda0, K)
    fun0 <- ll0 + joint.prior(beta, nu, Lw0, Aw0, U0, V0, lambda0, psi0)
    fun_seq <- c(fun0)
    ll_seq <- c(ll0)
  }
  beta_seq <- c(beta)
  time_seq <- c(0)
  start_time <- proc.time()[3]
  if (record_weights)
    w_seq <- list(Matrix::Matrix(w0, sparse = TRUE))
  if (verbose)
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
      ll <- joint.negloglikelihood(Lw, lambda, K)
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
    if (verbose)
      pb$tick(token = list(beta = beta, kth_eigval = eigvals[k], relerr = 2*max(werr/(w + w0), na.rm = 'ignore')))
    if (!fix_beta) {
      n_zero_eigenvalues <- sum(abs(eigvals) < eigtol)
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
  results <- list(laplacian = Lw, adjacency = Aw, w = w, psi = psi,
                  lambda = lambda, V = V, U = U, elapsed_time = time_seq,
                  beta_seq = beta_seq, convergence = has_w_converged)
  if (record_objective) {
    results$obj_fun <- fun_seq
    results$negloglike <- ll_seq
  }
  if (record_weights)
    results$w_seq <- w_seq
  return(results)
}
