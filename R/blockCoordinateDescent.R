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
  grad_f <- Lstar(L(w) - crossprod(sqrt(lambda) * t(U)) + Km / beta)
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
  # Updates the value of lambda
  #
  # Args:
  #   lb, ub: lower and upper bounds on the Lambda vector components
  #   β: scalar that controls the strength of the regularization term
  #   U: matrix
  #   w: vector
  #   N: dimension of each data sample
  #   K: number of components of the graph
  #
  # Returns:
  #   lambda: the updated value of lambda
  q <- N - K
  d <- diag(t(U) %*% L(w) %*% U)

  lambda <- .5 * (d + sqrt(d^2 + 4 / beta)) # unconstrained solution as initial point
  condition <- all(lambda[q] <= ub, lambda[1] >= lb,
                   lambda[2:q] >= lambda[1:(q-1)])

  while (!condition) {
    geq <- c(lb >= lambda[1], lambda[1:(q-1)] >= lambda[2:q], lambda[q] >= ub)
    flag1 <- geq[1]
    flag2 <- geq[q + 1]
    c1 <- 1
    c2 <- q + 1

    for (i in 1:q) {
      if (flag1) {
        lambda[i] <- lb
        flag1 <- geq[i + 1]
        if(!flag1) {
          c1 <- i + 1
        }
      }

      if (flag2) {
        lambda[q - i + 1] <- ub
        flag2 <- geq[q + 1 - i]
        if (!flag2) {
          c2 <- q + 1 - i
        }
      }

      if ((!flag1) & (!flag2)) {
        break
      }
    }

    if(c2 > c1) {
      if (c2 < (q + 1)) {
        geq3 <- lambda[c1:(c2-1)] >= lambda[(c1+1):c2]
      } else {
        geq3 <- lambda[c1:(c2-2)] >= lambda[(c1+1):(c2-1)]
      }
    } else {
      break
    }

    m <- c()
    for (i in 1:length(geq3)) {
        if (geq3[i]) {
          m <- c(m, c1 + i - 1)
        } else {
          d_mean <- mean(d[m])
          lambda[m] <- d_mean + sqrt(d_mean^2 + 4/beta)
          m <- c()
        }
    }
    condition <- all(lambda[q] <= ub, lambda[1] >= lb,
                     lambda[2:q] >= lambda[1:(q-1)])
  }

  return(lambda)
}
