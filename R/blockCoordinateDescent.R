# This module defines the updates of the block coordinate descent algorithm.
# More precisely, there are three variables to be updated:
#  1) a vector w (see expression (9)),
#  2) a matrix U (see expression (11)),
#  3) a matrix Λ (see expression (15)).

w_update <- function(w, U, beta, lambda, N, Km) {
  # Updates the value of w
  #
  # Args:
  #   w: vector
  #   U: matrix
  #   β: scalar
  #   Λ: vector
  #   n: dimension of each data sample
  #   Km: matrix (see section 1.1 for its definition)
  #
  # Returns:
  #   w_update: the updated value of w
  #grad_f <- Lstar(L(w) - U %*% diag(lambda) %*% t(U) + Km / beta)
  grad_f <- Lstar(L(w) - crossprod(t(U)*sqrt(lambda)) + Km / beta)
  w_update <- w - .5 * grad_f / N
  return(pmax(0, w_update))
}


U_update <- function(w, N, K) {
  # Updates the value of U
  #
  # Args:
  #   w: vector
  #   N: dimension of each data sample
  #   K: number of components
  #
  # Returns:
  #   U_update: the updated value of U
  return(eigen(L(w))$vectors[, (N-K):1])
}


lambda_update <- function(lb, ub, beta, U, w, N, K) {
  l <- N - K
  d <- diag(t(U) %*% L(w) %*% U)
  lambda <- .5 * (d + sqrt(d^2 + 4 / beta))
  condition <- all(lambda[l] <= ub, lambda[1] >= lb, lambda[2:l] >= lambda[1:(l-1)])
  if (condition) {
    return(lambda)
  }

  while (!condition) {
    geq <- c(lb >= lambda[1], lambda[1:(l-1)] >= lambda[2:l], lambda[l] >= ub)
    flag1 <- geq[1]
    flag2 <- geq[l]
    for (i in 1:l) {
      if (flag1) {
        lambda[i] <- lb
        flag1 <- geq[i+1]
      } else {
        c1 <- i
        flag1 <- FALSE
      }

      if (flag2) {
        lambda[l-i+1] <- ub
        flag2 <- geq[l-i]
      } else {
        c2 <- l-i+1
        flag2 <- FALSE
      }

      if ((!flag1) & (!flag2)) {
        break
      }
    }
    condition <- all(lambda[l] <= ub, lambda[1] >= lb, lambda[2:l] >= lambda[1:(l-1)])
  }
}
