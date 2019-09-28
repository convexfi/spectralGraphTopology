# implement the estimation algorithm proposed by Elgimez and Ortega 2017
# named Combinatorial Graph Laplacian

learn_combinatorial_graph_laplacian <- function(S, A_mask, alpha, prob_tol = 1e-4,
                                                inner_tol = 1e-4, maxiter = 1000, regtype = 1) {
  n <- nrow(S)
  e_v <- rep(1, n) / sqrt(n)
  dc_var <- t(e_v) %*% S %*% e_v
  isshifting <- c(abs(dc_var) < prob_tol)
    
  if (isshifting) {
      S <- S + 1 / n
  }
  if (regtype == 1) { 
      H <- alpha * (2 * diag(n) - matrix(1, n, n))
  } else if (regtype == 2) {
      H <- alpha * (diag(n) - matrix(1, n, n))
  }
  K <- S + H
  O_init <- diag(1 / diag(K))
  C <- diag(diag(K))
  O <- O_init

  O_best <- O
  C_best <- C

  frob_norm <- c()
  time_counter <- c()
  has_converged <- FALSE
  cycle <- 0

  for (i in c(1:maxiter)) {
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

    frob_norm <- c(frob_norm, norm(O_old - O, 'F') / norm(O_old, "F"))
    if (i > 6) {
      if (frob_norm[i] < prob_tol) {
        O_best <- O
        C_best <- C
        break
      }
    }
  }
  O <- O_best - (1 / n)
  C <- C_best - (1 / n)
  return(list(Theta = O, Cov = C, frob_norm = frob_norm))
}

update_sherman_morrison_diag <- function(O, C, shift, idx) {
  O[idx, idx] <- O[idx, idx] + shift
  c_d <- C[idx, idx]
  C <- C - (shift / (1 + shift * c_d)) * C[, idx] %*% t(C[idx, ])
  return(list(O = O, C = C))
}



