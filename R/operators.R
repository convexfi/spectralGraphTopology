# This module defines two operators, namely, L and Lstar,
# used in the graph learning problem.
# I'm trying to follow Google R's style guide:
# https://google.github.io/styleguide/Rguide.xml


LOp <- function(w, n) {
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
  Lw <- matrix(0, n, n)
  k <- length(w)
  for(i in ((n-1):1)) {
    j <- n - i
    Lw[i, (i+1):n] = -tail(w[1:k], j)
    k <- k - j
  }

  Lw <- Lw + t(Lw)
  Lw <- Lw - diag(colSums(Lw))

  return(Lw)
}


LStarOp <- function(W) {
  # The L star operator, LStarOp, acts upon a matrix W and
  # returns a vector given by expression (7) in the Work Plan I.
  #
  # Args:
  #   W: matrix which LStarOp will act upon
  #
  # Returns:
  #   w: vector

  return(NA)
}
