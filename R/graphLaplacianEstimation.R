# this module contains implementations of the algorithms
# proposed in L. Zhao et. al. "Optimization Algorithms for
# Graph Laplacian Estimation via ADMM and MM", IEEE Trans. Sign. Proc. 2019.

get_incidence_from_adjacency <- function(A) {
  p <- nrow(A)
  m <- .5 * sum(A > 0)
  E <- matrix(0, p, m)
  k <- 1
  for (i in c(1:(p-1))) {
    for (j in c((i+1):p)) {
      if (A[i, j] > 0) {
        E[i, k] <- 1
        E[j, k] <- -1
        k <- k + 1
      }
    }
  }
  return(E)
}

#' @title Learn the weighted Laplacian matrix of a graph using the MM method
#'
#' @param S a pxp sample covariance/correlation matrix
#' @param A the binary adjacency matrix of the graph
#' @param alpha L1 regularization hyperparameter
#' @param maxiter the maximum number of iterations
#' @param reltol relative tolerance on the weight vector w
#' @param abstol absolute tolerance on the weight vector w
#' @param record_objective whether or not to record the objective function. Default is FALSE
#' @param verbose if TRUE, then a progress bar will be displayed in the console. Default is TRUE
#' @return A list containing possibly the following elements:
#' \item{\code{Laplacian}}{the estimated Laplacian Matrix}
#' \item{\code{Adjacency}}{the estimated Adjacency Matrix}
#' \item{\code{convergence}}{boolean flag to indicate whether or not the optimization converged}
#' \item{\code{obj_fun}}{values of the objective function at every iteration in case record_objective = TRUE}
#' @author Ze Vinicius, Jiaxi Ying, and Daniel Palomar
#' @references Licheng Zhao, Yiwei Wang, Sandeep Kumar, and Daniel P. Palomar.
#'             Optimization Algorithms for Graph Laplacian Estimation via ADMM and MM.
#'             IEEE Trans. on Signal Processing, vol. 67, no. 16, pp. 4231-4244, Aug. 2019

#' @export
learn_laplacian_gle_mm <- function(S, A, alpha = 0, maxiter = 10000, reltol = 1e-5,
                                   abstol = 1e-5, record_objective = FALSE,
                                   verbose = TRUE) {
  Sinv <- MASS::ginv(S)
  mask <- Ainv(A) > 0
  w <- w_init("naive", Sinv)[mask]
  wk <- w
  # number of nodes
  p <- nrow(S)
  # number of nonzero edges
  m <- .5 * sum(A > 0)
  # l1-norm penalty factor
  J <- matrix(1, p, p) / p
  H <- 2 * diag(p) - p * J
  K <- S + alpha * H
  E <- get_incidence_from_adjacency(A)
  R <- t(E) %*% K %*% E
  r <- nrow(R)
  G <- cbind(E, rep(1, p))
  if (record_objective)
    fun <- obj_func(E, K, wk, J)
  if (verbose)
    pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta",
                                     total = maxiter, clear = FALSE, width = 80)
  for (k in c(1:maxiter)) {
    w_aug <- c(wk, 1 / p)
    G_aug_t <- t(G) * w_aug
    G_aug <- t(G_aug_t)
    Q <- G_aug_t %*% solve(G_aug %*% t(G), G_aug)
    Q <- Q[1:m, 1:m]
    wk <- sqrt(diag(Q) / diag(R))
    if (record_objective)
      fun <- c(fun, obj_func(E, K, wk, J))
    if (verbose)
       pb$tick()
    werr <- abs(w - wk)
    has_converged <- all(werr <= .5 * reltol * (w + wk)) || all(werr <= abstol)
    if (has_converged && k > 1) break
    w <- wk
  }
  z <- rep(0, .5 * p * (p - 1))
  z[mask] <- wk
  results <- list(Laplacian = L(z), Adjacency = A(z), maxiter = k,
                  convergence = has_converged)
  if (record_objective)
    results$obj_fun <- fun
  return(results)
}

obj_func <- function(E, K, w, J) {
  p <- ncol(J)
  EWEt <- E %*% diag(w) %*% t(E)
  Gamma <- EWEt + J
  lambda <- eigval_sym(Gamma)[2:p]
  return(sum(diag(E %*% diag(w) %*% t(E) %*% K)) - sum(log(lambda)))
}

#' @title Learn the weighted Laplacian matrix of a graph using the ADMM method
#'
#' @param S a pxp sample covariance/correlation matrix
#' @param A the binary adjacency matrix of the graph
#' @param alpha L1 regularization hyperparameter
#' @param rho ADMM convergence rate hyperparameter
#' @param maxiter the maximum number of iterations
#' @param reltol relative tolerance on the Laplacian matrix estimation
#' @param record_objective whether or not to record the objective function. Default is FALSE
#' @param verbose if TRUE, then a progress bar will be displayed in the console. Default is TRUE
#' @return A list containing possibly the following elements:
#' \item{\code{Laplacian}}{the estimated Laplacian Matrix}
#' \item{\code{Adjacency}}{the estimated Adjacency Matrix}
#' \item{\code{convergence}}{boolean flag to indicate whether or not the optimization converged}
#' \item{\code{obj_fun}}{values of the objective function at every iteration in case record_objective = TRUE}
#' @author Ze Vinicius, Jiaxi Ying, and Daniel Palomar
#' @references Licheng Zhao, Yiwei Wang, Sandeep Kumar, and Daniel P. Palomar.
#'             Optimization Algorithms for Graph Laplacian Estimation via ADMM and MM.
#'             IEEE Trans. on Signal Processing, vol. 67, no. 16, pp. 4231-4244, Aug. 2019
#' @export
learn_laplacian_gle_admm <- function(S, A, alpha = 0, rho = 1, maxiter = 10000,
                                     reltol = 1e-5, record_objective = FALSE,
                                     verbose = TRUE) {
  p <- nrow(S)
  Sinv <- MASS::ginv(S)
  w <- w_init("naive", Sinv)
  Yk <- L(w)
  Theta <- Yk
  P <- eigvec_sym(Yk)
  P <- P[, 2:p]
  Ck <- Yk
  # l1-norm penalty factor
  J <- matrix(1, p, p) / p
  H <- 2 * diag(p) - p * J
  K <- S + alpha * H
  if (verbose)
    pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta",
                                     total = maxiter, clear = FALSE, width = 80)
  if (record_objective)
    fun <- c()
  # ADMM loop
  for (k in c(1:maxiter)) {
    Gamma <- t(P) %*% (K + Yk - rho * Ck) %*% P
    U <- eigvec_sym(Gamma)
    lambda <- eigval_sym(Gamma)
    d <- .5 * c(sqrt(lambda ^ 2 + 4 / rho) - lambda)
    Xik <- crossprod(sqrt(d) * t(U))
    Thetak <- P %*% Xik %*% t(P)
    Ck <- (diag(pmax(0, diag(Yk / rho + Thetak))) +
           A * pmin(0, Yk / rho + Thetak))
    Yk <- Yk + rho * (Thetak - Ck)
    if (record_objective)
      fun <- c(fun, aug_lag(K, P, Xik, Yk, Ck, d, rho))
    has_converged <- norm(Theta - Thetak, "F") / norm(Thetak, "F") < reltol
    if (has_converged && k > 1) break
    Theta <- Thetak
    if (verbose)
      pb$tick()
  }
  results <- list(Laplacian = Thetak, Adjacency = diag(diag(Thetak)) - Thetak,
                  convergence = has_converged)
  if (record_objective)
    results$obj_fun <- fun
  return(results)
}

# compute the partial augmented Lagragian
aug_lag <- function(K, P, Xi, Y, C, d, rho) {
  PXiPt <- P %*% Xi %*% t(P)
  return(sum(diag(Xi %*% t(P) %*% K %*% P)) - sum(log(d)) +
         sum(diag(t(Y) %*% (PXiPt - C))) + .5 * rho * norm(PXiPt - C, "F") ^ 2)
}
