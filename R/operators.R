# This module defines two operators, namely, L and Lstar,
# used in the graph learning problem.
# I'm trying to follow Google R's style guide:
# https://google.github.io/styleguide/Rguide.xml


LOp <- function(w) {
  # The L operator, LOp, acts upon a vector w and returns
  # a matrix given by expression (5) in the Work Plan I.
  #
  # Args:
  #   w: vector which LOp will act upon.
  #   n: dimension of one data sample, which is the order
  #      of the matrix Lw
  #
  # Returns:
  #   Lw: matrix
  k <- length(w)
  n <- as.integer(.5 * (1 + sqrt(1 + 8 * k)))
  Lw <- matrix(0, n, n)
  for (i in ((n-1):1)) {
    j <- n - i
    Lw[i, (i+1):n] <- -tail(w[1:k], j)
    k <- k - j
  }

  Lw <- Lw + t(Lw)
  Lw <- Lw - diag(colSums(Lw))

  return(Lw)
}

LStarOp <- function(Y) {
  # The L star operator, LStarOp, acts upon a matrix W and
  # returns a vector given by expression (7) in the Work Plan I.
  #
  # Args:
  #   Y: matrix which LStarOp will act upon
  #
  # Returns:
  #   LStarY: vector
  N <- ncol(Y)
  k <- as.integer(N * (N - 1) / 2)
  LStarY <- array(0., k)
  j <- 1
  l <- 2
  for (i in 1:k) {
    LStarY[i] <- Y[j, j] + Y[l, l] - Y[l, j] - Y[j, l]
    if (l == N) {
      j <- j + 1
      l <- j + 1
    } else {
      l <- l + 1
    }
  }
  return(LStarY)
}


LStarOpImpl <- function(Y) {
  # An alternative implementation of the L star operator.
  # Implemented only for the purposes of testing.
  #
  # Args:
  #   Y: matrix which LStarOp will act upon
  #
  # Returns:
  #   LStarY: vector
  n <- ncol(Y)
  k <- as.integer(n * (n - 1) / 2)
  LStarY <- array(0., k)
  for (i in 1:k) {
    w <- array(0., k)
    w[i] <- 1.
    Lw <- LOp(w)
    LStarY[i] <- sum(diag(t(Y) %*% Lw))
  }

  return(LStarY)
}
