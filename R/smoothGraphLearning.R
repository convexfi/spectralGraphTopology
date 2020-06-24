# this module contains implementations of graph learning methods
# from smooth signals

#' @title Learn graphs from a smooth signal representation approach
#'
#' This function learns a graph from a observed data matrix using the
#' method proposed by Dong (2016).
#'
#' @param X a p-by-n data matrix, where p is the number of nodes and n is the
#'        number of observations
#' @param alpha hyperparameter that controls the importance of the Dirichlet
#'        energy penalty
#' @param beta hyperparameter that controls the importance of the L2-norm
#'        regularization
#' @param maxiter maximum number of iterations
#' @param ftol relative error on the objective function to be used as the
#'        stopping criteria
#' @param verbose if TRUE, then a progress bar will be displayed in the console. Default is TRUE
#' @return A list containing the following items
#' \item{\code{laplacian}}{estimated Laplacian Matrix}
#' \item{\code{Y}}{a smoothed approximation of the data matrix X}
#' \item{\code{convergence}}{whether or not the algorithm has converged within the tolerance and max number of iterations}
#' \item{\code{obj_fun}}{objective function value at every iteration, in case record_objective = TRUE}
#' @references X. Dong, D. Thanou, P. Frossard and P. Vandergheynst, "Learning
#'             Laplacian Matrix in Smooth Graph Signal Representations,"
#'             in IEEE Transactions on Signal Processing, vol. 64, no. 23,
#'             pp. 6160-6173, Dec.1, 2016.
#' @export
learn_graph_sigrep <- function(X, alpha = 1e-3, beta = 1e-1, maxiter = 1000, ftol = 1e-4, verbose = TRUE) {
  p <- nrow(X)
  Y <- X
  obj_values <- c()
  fun0 <- Inf
  if (verbose)
    pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta  relerr: :relerr",
                                     total = maxiter, clear = FALSE, width = 80)
  for (k in c(1:maxiter)) {
    L <- glsigrep.update_L(Y, alpha, beta, p)
    Y <- glsigrep.update_Y(L, X, alpha, p)
    funk <- glsigrep.obj_function(X, Y, L, alpha, beta)
    obj_values <- c(obj_values, funk)
    relerr <- abs(funk - fun0) / fun0
    if (k > 1 && abs(funk - fun0) / fun0 < ftol) break
    if (verbose) pb$tick(token = list(relerr = relerr))
    fun0 <- funk
  }
  return(list(laplacian = L, Y = Y, convergence = (k < maxiter), obj_fun = obj_values))
}


glsigrep.update_L <- function(Y, alpha, beta, p) {
  L <- CVXR::Symmetric(p)
  ones <- rep(1, p)
  zeros <- rep(0, p)
  zeros_mat <- matrix(0, p, p)
  objective <- CVXR::Minimize(alpha * CVXR::matrix_trace(L %*% Y %*% t(Y))
                              + beta * CVXR::sum_squares(L))
  constraints <- list(CVXR::matrix_trace(L) == p,
                      L %*% ones == zeros,
                      CVXR::upper_tri(L) <= 0,
                      CVXR::diag(L) > 0)
  problem <- CVXR::Problem(objective, constraints)
  result <- solve(problem)
  return(as.matrix(result$getValue(L)))
}

glsigrep.update_Y <- function(L, X, alpha, p) {
  Ip <- diag(p)
  return(solve(Ip + alpha * L, X))
}

glsigrep.obj_function <- function(X, Y, L, alpha, beta) {
  return(norm(X - Y, 'F') ^ 2 + alpha * sum(diag(L %*% Y %*% t(Y))) + beta * norm(L, 'F') ^ 2)
}


#' @title Learn a graph from smooth signals
#'
#' This function learns a connected graph given an observed signal matrix
#' using the method proposed by Kalofilias (2016).
#'
#' @param X a p-by-n data matrix, where p is the number of nodes and n is the
#'        number of observations
#' @param alpha hyperparameter that controls the importance of the Dirichlet
#'        energy penalty
#' @param beta hyperparameter that controls the importance of the L2-norm
#'        regularization
#' @references V. Kalofolias, "How to learn a graph from smooth signals", in Proc. Int.
#'             Conf. Artif. Intell. Statist., 2016, pp. 920–929.
#' @export
learn_smooth_graph <- function(X, alpha = 1e-2, beta = 1e-4, step_size = 1e-2,
                               maxiter = 1000, tol = 1e-4) {
  p <- nrow(X)
  S <- Sop(p)
  wk <- spectralGraphTopology:::w_init("naive", MASS::ginv(cor(t(X))))
  dk <- D(wk)
  ## constants
  mu <- 2 * beta + sqrt(2 * (p - 1))
  eps <- lin_map(0, 0, 1 / (1 + mu))
  gamma <- lin_map(step_size, eps, (1 - eps) / mu)
  z <- upper_view_vec(pairwise_matrix_rownorm2(X))
  for (k in 1:maxiter) {
    wk_prev <- wk
    dk_prev <- dk
    yk <- wk - gamma * (2 * beta * wk + t(S) %*% dk)
    yk_bar <- dk + gamma * S %*% wk
    pk <- pmax(0, yk - 2 * gamma * z)
    pk_bar <- .5 * (yk_bar - sqrt(yk_bar * yk_bar + 4 * alpha * gamma))
    qk <- pk - gamma * (2 * beta * pk + t(S) %*% pk_bar)
    qk_bar <- pk_bar + gamma * S %*% pk
    wk <- wk - yk + qk
    dk <- dk - yk_bar + qk_bar
    if(norm(wk - wk_prev, '2') / norm(wk_prev, '2') < tol &&
       norm(dk - dk_prev, '2') / norm(dk_prev, '2') < tol)
      break
  }
  return(list(laplacian = L(wk), adjacency = A(wk), convergence = (k < maxiter)))
}

Sop <- function(p) {
  ncols <- .5 * p * (p - 1)
  ii <- rep(0, ncols)
  jj <- rep(0, ncols)
  S <- matrix(0, p, ncols)

  k <- 1
  for (i in c(2:p)) {
    ii[k:(k + (p - i))] <- c(i:p)
    k <- k + (p - i + 1)
  }
  k <- 1
  for (i in c(2:p)) {
    jj[k:(k + p - i)] <- i - 1
    k <- k + (p - i + 1)
  }
  r <- c(c(1:ncols), c(1:ncols))
  c <- c(ii, jj)
  for (k in c(1:length(r)))
    S[c[k], r[k]] <- 1
  return(S)
}

lin_map <- function(x, a, b) {
  return(x * (b - a) + a)
}
