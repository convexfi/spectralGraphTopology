# This module defines the updates of the block coordinate descent algorithm.
# More precisely, there are three variables to be updated:
#  1) a vector w (see expression (9)),
#  2) a matrix U (see expression (11)),
#  3) a matrix Λ (see expression (15)).

w_update <- function(w, U, beta, Lambda, n, Km) {
  # Function to update the value of w
  #
  # Args:
  #   w: vector
  #   U: matrix
  #   β: scalar
  #   Λ: matrix
  #   n: dimension of each data sample
  #   Km: matrix (see section 1.1 for its definition)
  #
  # Returns:
  #   w_update: the updated value of w

  grad_f <- LStarOp(LOp(w)) - LStarOp(U %*% Lambda %*% t(U) - Km / beta)
  w_update <- w - .5 * grad_f / n
  mask <- w_update < 0
  w_update[mask] <- 0
  return(w_update)
}
