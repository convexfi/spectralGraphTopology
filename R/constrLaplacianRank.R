constr_laplacian_rank <- function(Y, k = 1, m = 5, S0 = NULL, lmd = 1, eig_tol = 1e-9,
                                  edge_tol = 1e-6, maxiter = 1000, regularization_type = 2) {
  time_seq <- c(0)
  start_time <- proc.time()[3]
  A <- build_initial_graph(Y, m)
  n <- ncol(A)
  if (is.null(S0))
    S <- matrix(1/n, n, n)
  else
    S <- S0
  DS <- diag(.5 * colSums(S + t(S)))
  LS <-  DS - .5 * (S + t(S))
  DA <- diag(.5 * colSums(A + t(A)))
  LA <- DA - .5 * (A + t(A))
  #if (k == 1)
  F <- matrix(eigvec_sym(LA)[, 1:k])
  #else
  #  F <- eigvec_sym(LA)[, 1:k]
  # bounds for variables in the QP solver
  bvec <- c(1, rep(0, n))
  Amat <- cbind(rep(1, n), diag(n))
  lmd_seq <- c(lmd)
  pb <- progress::progress_bar$new(format = "<:bar> :current/:total  eta: :eta  lambda: :lmd  null_eigvals: :null_eigvals",
                                   total = maxiter, clear = FALSE, width = 100)
  for (ii in c(1:maxiter)) {
    V <- pairwise_matrix_rownorm(F)
    if (regularization_type == 1) {
      for (i in c(1:n)) {
        u <- .5 / abs(S[i, ] - A[i, ])
        p <- u * A[i, ] - .5 * lmd * V[i, ]
        qp <- quadprog::solve.QP(Dmat = diag(u), dvec = p, Amat = Amat, bvec = bvec, meq = 1)
        S[i, ] <- qp$solution
      }
    } else if (regularization_type == 2) {
      for (i in c(1:n)) {
        p <- A[i, ] - .5 * lmd * V[i, ]
        qp <- quadprog::solve.QP(Dmat = diag(n), dvec = p, Amat = Amat, bvec = bvec, meq = 1)
        S[i, ] <- qp$solution
      }
    }
    DS <- diag(.5 * colSums(S + t(S)))
    LS <- DS - .5 * (S + t(S))
    F <- eigvec_sym(LS)[, 1:k]
    eig_vals <- eigval_sym(LS)
    n_zero_eigenvalues <- sum(abs(eig_vals) < eig_tol)
    time_seq <- c(time_seq, proc.time()[3] - start_time)
    pb$tick(token = list(lmd = lmd, null_eigvals = n_zero_eigenvalues))
    if (k < n_zero_eigenvalues)
      lmd <- .5 * lmd
    else if (k > n_zero_eigenvalues)
      lmd <- 2 * lmd
    else
      break
    lmd_seq <- c(lmd_seq, lmd)
  }
  LS[abs(LS) < edge_tol] <- 0
  AS <- diag(diag(LS)) - LS
  return(list(Laplacian = LS, Adjacency = AS, eigenvalues = eig_vals,
              lmd_seq = lmd_seq, elapsed_time = time_seq))
}


build_initial_graph <- function(Y, m) {
  n <- nrow(Y)
  A <- matrix(0, n, n)
  E <- pairwise_matrix_rownorm(Y)
  for (i in c(1:n)) {
    sorted_index <- order(E[i, ])
    j_sweep <- sorted_index[2:(m+1)]
    den <- m * E[i, sorted_index[m+2]] - sum(E[i, j_sweep])
    ei <- E[i, sorted_index[m+2]]
    for (j in j_sweep) {
      A[i, j] <- (ei - E[i, j]) / den
    }
  }
  return(A)
}


objective_function <- function(A, S, LS, F, lmd) {
  return(sum(abs(A - S)) + 2 * lmd * sum(diag(t(F) %*% LS %*% F)))
}