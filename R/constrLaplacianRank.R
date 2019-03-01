constr_laplacian_rank <- function(Y, k = 1, m = 5, lmd = 1e-1, eig_tol = 1e-4,
                                  maxiter = 1000, regularization_type = 1) {
  A <- initial_graph(Y, m)
  n <- ncol(A)
  S <- matrix(1, n, n) / n
  D <- diag(.5 * colSums(A + t(A)))
  F <- eigenvectors(D - .5 * (A + t(A)))[, 1:k]
  # bounds for variables in the QP solver
  bvec <- c(1, rep(0, n))
  Amat <- cbind(rep(1, n), diag(n))
  # objective function placeholder
  fun_seq <- c()
  fun_k <- Inf
  pb = txtProgressBar(min = 0, max = maxiter, initial = 0)
  for (ii in c(1:maxiter)) {
    setTxtProgressBar(pb, ii)
    if (regularization_type == 1) {
      V <- pairwise_matrix_rownorm(F)
      for (i in c(1:n)) {
        U <- diag(.5 / abs(S[i, ] - A[i, ]))
        p <- U %*% A[i, ] - .5 * lmd * V[i, ]
        qp <- quadprog::solve.QP(Dmat = U, dvec = p, Amat = Amat, bvec = bvec, meq = 1)
        S[i, ] <- qp$solution
      }
    } else if (regularization_type == 2) {

    } else {
      stop("type ", type, " not recognized.")
    }
    D <- diag(.5 * colSums(S + t(S)))
    LS <- D - .5 * (S + t(S))
    if (k < sum(abs(eigenvalues(LS)) < .Machine$double.eps))
      lmd <- .5 * lmd
    else
      lmd <- 2 * lmd
    F <- eigenvectors(LS)[, 1:k]
    criteria <- sum(abs(eigenvalues(S)[1:k]))
    print(criteria)
    fun_next <- objective_function(A, S, LS, F, lmd)
    fun_seq <- c(fun_seq, fun_next)
    if (criteria < k * eig_tol)
      break
  }
  print(fun_seq)
  return(LS)
}

objective_function <- function(A, S, LS, F, lmd) {
  return(sum(abs(A - S)) + 2 * lmd * sum(diag(t(F) %*% LS %*% F)))
}


initial_graph <- function(Y, m) {
  n <- nrow(Y)
  A <- matrix(0, n, n)
  E <- pairwise_matrix_rownorm(Y)
  for (i in c(1:n)) {
    sorted_index <- order(E[i, ])
    for (j in sorted_index[2:m]) {
      A[i, j] <- (E[i, sorted_index[m+1]] - E[i, j]) / (m * E[i, sorted_index[m+1]] - sum(E[i, sorted_index[2:m]]))
      A[j, i] <- A[i, j]
    }
  }
  return(A)
}


# given a n-by-m matrix M, this function outputs a n-by-n matrix V such
# that V[i, j] = ||M[i, ] - M[j, ]||_2 ^ 2
pairwise_matrix_rownorm <- function(M) {
  n <- nrow(M)
  V <- matrix(0, n, n)
  for (i in c(1:(n-1)))
    for (j in c((i+1):n)) {
      V[i, j] = norm((M[i, ] - M[j, ]), "2") ^ 2
      V[j, i] <- V[i, j]
    }
  return(V)
}
