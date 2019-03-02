constr_laplacian_rank <- function(Y, k = 1, m = 5, S0 = NULL, lmd = 1, eig_tol = 1e-6,
                                  edge_tol = 1e-6, maxiter = 1000, regularization_type = 2) {
  A <- build_initial_graph(Y, m)
  n <- ncol(A)
  if (is.null(S0)) S <- matrix(1/n, n, n)
  else S <- S0
  DS <- diag(.5 * colSums(S + t(S)))
  LS <-  DS - .5 * (S + t(S))
  DA <- diag(.5 * colSums(A + t(A)))
  LA <- DA - .5 * (A + t(A))
  F <- eigenvectors(LA)[, 1:k]
  # bounds for variables in the QP solver
  bvec <- c(1, rep(0, n))
  Amat <- cbind(rep(1, n), diag(n))
  # objective function placeholder
  fun_k <- objective_function(A, S, LS, F, lmd)
  fun_seq <- c(fun_k)
  lmd_seq <- c(lmd)
  pb = txtProgressBar(min = 0, max = maxiter, initial = 0)
  for (ii in c(1:maxiter)) {
    setTxtProgressBar(pb, ii)
    V <- pairwise_matrix_rownorm(F)
    if (regularization_type == 1) {
      for (i in c(1:n)) {
        U <- diag(.5 / abs(S[i, ] - A[i, ]))
        p <- diag(U) * A[i, ] - .5 * lmd * V[i, ]
        qp <- quadprog::solve.QP(Dmat = U, dvec = p, Amat = Amat, bvec = bvec, meq = 1)
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
    F <- eigenvectors(LS)[, 1:k]
    fun_next <- objective_function(A, S, LS, F, lmd)
    fun_seq <- c(fun_seq, fun_next)
    eig_vals <- eigenvalues(LS)
    n_zero_eigenvalues <- sum(abs(eig_vals) < eig_tol)
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
              obj_fun = fun_seq, lmd_seq = lmd_seq))
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

