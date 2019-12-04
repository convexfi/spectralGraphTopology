# implement the estimation algorithm proposed by Elgimez and Ortega 2017
# named Combinatorial Graph Laplacian

#' @title Learn the Combinatorial Graph Laplacian from data
#'
#' Learns a graph Laplacian matrix using the Combinatorial Graph Laplacian (CGL)
#' algorithm proposed by Egilmez et. al. (2017)
#'
#' @param S sample covariance matrix
#' @param A_mask binary adjacency matrix of the graph
#' @param alpha L1-norm regularization hyperparameter
#' @param reltol minimum relative error considered for the stopping criteri
#' @param max_cycle maximum number of cycles
#' @param regtype type of L1-norm regularization. If reg_type == 1, then all
#'        elements of the Laplacian matrix will be regularized. If reg_type == 2,
#'        only the off-diagonal elements will be regularized
#' @param record_objective whether or not to record the objective function value
#'        at every iteration. Default is FALSE
#' @param verbose if TRUE, then a progress bar will be displayed in the console. Default is TRUE
#' @return A list containing possibly the following elements
#' \item{\code{Laplacian}}{estimated Laplacian Matrix}
#' \item{\code{elapsed_time}}{elapsed time recorded at every iteration}
#' \item{\code{frod_norm}}{relative Frobenius norm between consecutive estimates of the Laplacian matrix}
#' \item{\code{convergence}}{whether or not the algorithm has converged within the tolerance and max number of iterations}
#' \item{\code{obj_fun}}{objective function value at every iteration, in case record_objective = TRUE}
#' @references H. E. Egilmez, E. Pavez and A. Ortega, "Graph Learning From Data
#'             Under Laplacian and Structural Constraints", in IEEE Journal of
#'             Selected Topics in Signal Processing, vol. 11, no. 6, pp. 825-841, Sept. 2017.
#'             Original MATLAB source code is available at: https://github.com/STAC-USC/Graph_Learning
#' @export
learn_combinatorial_graph_laplacian <- function(S, A_mask = NULL, alpha = 0, reltol = 1e-5,
                                                max_cycle = 10000, regtype = 1,
                                                record_objective = FALSE, verbose = TRUE) {
  n <- nrow(S)
  if (is.null(A_mask))
    A_mask <- matrix(1, n, n) - diag(n)
  e_v <- rep(1, n) / sqrt(n)
  dc_var <- t(e_v) %*% S %*% e_v
  isshifting <- c(abs(dc_var) < reltol)
  if (isshifting) {
    S <- S + 1 / n
  }
  if (regtype == 1) {
    H <- alpha * (diag(n) - matrix(1, n, n))
  } else if (regtype == 2) {
    H <- alpha * (2 * diag(n) - matrix(1, n, n))
  }
  K <- S + H
  O_init <- diag(1 / diag(K))
  C <- diag(diag(K))
  O <- O_init
  O_best <- O
  C_best <- C
  frob_norm <- c()
  has_converged <- FALSE
  if (verbose)
    pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta",
                                     total = max_cycle, clear = FALSE, width = 80)
  time_seq <- c(0)
  if (record_objective)
    fun <- vanilla.objective(O_best - (1/n), K)
  start_time <- proc.time()[3]
  for (i in c(1:max_cycle)) {
    O_old <- O
    for (u in c(1:n)) {
      minus_u <- setdiff(c(1:n), u)
      k_u <- K[minus_u, u]
      k_uu <- K[u, u]
      c_u <- C[minus_u, u]
      c_uu <- C[u, u]
      Ou_i <- C[minus_u, minus_u] - (c_u %*% t(c_u) / c_uu)
      # block-descent variables
      beta <- rep(0, n-1)
      ind_nz <- A_mask[minus_u, u] == 1
      A_large <- Ou_i
      A_nnls <- Ou_i[ind_nz, ind_nz]
      b <- k_u / k_uu + (A_large %*% rep(1, n-1) / n)
      b_nnls <- b[ind_nz]
      # block-descent step
      Dmat <- A_nnls
      dvec <- b_nnls
      Amat <- diag(length(dvec))
      bvec <- rep(0, length(dvec))
      beta_quad <- - quadprog::solve.QP(Dmat = Dmat, dvec = dvec,
                                        Amat = Amat, bvec = bvec)$solution
      beta[ind_nz] <- beta_quad
      o_u <- beta + 1/n
      o_uu <- 1/k_uu + t(o_u) %*% Ou_i %*% o_u
      # Update the current Theta
      O[u, u] <- o_uu
      O[minus_u, u] <- o_u
      O[u, minus_u] <- o_u
      # Update the current Theta inverse
      cu <- (Ou_i %*% o_u) / c(o_uu - t(o_u) %*% Ou_i %*% o_u)
      cuu <- 1 / c(o_uu - t(o_u) %*% Ou_i %*% o_u)
      C[u, u] <- cuu
      C[u, minus_u] <- -cu
      C[minus_u, u] <- -cu
      C[minus_u, minus_u] <- (Ou_i + (cu %*% t(cu)) / cuu)
    }
    if (i > 4) {
      d_shifts <- O %*% rep(1, n) - 1
      large_diag_idx <- c(1:n)[abs(d_shifts) > 1e-12]
      for (idx_t in 1:length(large_diag_idx)) {
        idx <- large_diag_idx[idx_t]
        smd <- update_sherman_morrison_diag(O, C, -d_shifts[idx], idx)
        O <- smd$O
        C <- smd$C
      }
    }
    O_best <- O
    C_best <- C
    if (record_objective)
      fun <- c(fun, vanilla.objective(O_best - (1/n), K))
    if (verbose)
      pb$tick()
    time_seq <- c(time_seq, proc.time()[3] - start_time)
    frob_norm <- c(frob_norm, norm(O_old - O, 'F') / norm(O_old, "F"))
    if (i > 6) {
      if (frob_norm[i] < reltol) {
        O_best <- O
        C_best <- C
        break
      }
    }
  }
  if (i < max_cycle) has_converged = TRUE
  else has_converged = FALSE
  O <- O_best - (1 / n)
  C <- C_best - (1 / n)
  results <- list(Laplacian = O, frob_norm = frob_norm,
                  elapsed_time = time_seq, convergence = has_converged)
  if (record_objective)
    results$obj_fun <- fun
  return(results)
}

update_sherman_morrison_diag <- function(O, C, shift, idx) {
  O[idx, idx] <- O[idx, idx] + shift
  c_d <- C[idx, idx]
  C <- C - (shift / (1 + shift * c_d)) * C[, idx] %*% t(C[idx, ])
  return(list(O = O, C = C))
}
