constr_laplacian_rank <- function(A, k = 1, lmd = 1e3, eig_tol = 1e-4,
                                  maxiter = 1000, regularization_type = 1) {
  n <- ncol(A)
  S <- matrix(0, n, n)
  D <- diag(.5 * colSums(A + t(A)))
  F <- eigenvectors(D - .5 * (A + t(A)))[, 1:k]
  # bounds for variables in the QP solver
  bvec <- c(1, rep(0, n))
  Amat <- cbind(rep(1, n), diag(n))

  pb = txtProgressBar(min = 0, max = maxiter, initial = 0)
  for (ii in c(1:maxiter)) {
    setTxtProgressBar(pb, ii)
    if (regularization_type == 1) {
      V <- pairwise_matrix_row_norm(F)
      for (i in c(1:n)) {
        U <- diag(.5 / abs(S[i, ] - A[i, ]))
        p <- U %*% A[i, ] - .5 * lmd * V[i, ]
        #qp <- osqp::solve_osqp(U, p, Amat, l, u)
        #                       #osqp::osqpSettings(verbose = FALSE))
        qp <- quadprog::solve.QP(Dmat = U, dvec = p, Amat = Amat, bvec = bvec, meq = 1)
        S[i, ] <- qp$solution
      }
    } else if (regularization_type == 2) {
      stop("type ", type, " not implemented yet.")
    } else {
      stop("type ", type, " not recognized.")
    }
    D <- diag(.5 * colSums(S + t(S)))
    F <- eigenvectors(D - .5 * (S + t(S)))[, 1:k]
    criteria <- sum(abs(eigenvalues(S)[1:k]))
    print(criteria)
    if (criteria < k * eig_tol)
      break
  }
  return(S)
}

pairwise_matrix_row_norm <- function(M) {
  n <- nrow(M)
  V <- matrix(0, n, n)
  for (i in c(1:n))
    for (j in c(1:n))
      V[i, j] = norm((M[i, ] - M[j, ]), "2") ^ 2
  return(V + t(V))
}
