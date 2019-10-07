#' Constructs a block diagonal matrix from a list of square matrices
#'
#' @param ... list of matrices or individual matrices
#' @return block diagonal matrix
#' @examples
#' library(spectralGraphTopology)
#' X <- L(c(1, 0, 1))
#' Y <- L(c(1, 0, 1, 0, 0, 1))
#' B <- block_diag(X, Y)
#' B
#' @export
block_diag <- function(...) {
  return(blockDiagCpp(list(...)))
}


#' Computes the relative error between two matrices
#'
#' @param A first matrix
#' @param B second matrix
#' @examples
#' library(spectralGraphTopology)
#' X <- L(c(1, 0, 1))
#' relative_error(X, X)
#' @export
relative_error <- function(A, B) {
  return(norm(A - B, type = "F") / norm(B, type = "F"))
}


#' Computes the fscore between two matrices
#'
#' @param A first matrix
#' @param B second matrix
#' @param eps real number such that edges whose values are smaller than eps are
#'            not considered in the computation of the fscore
#' @examples
#' library(spectralGraphTopology)
#' X <- L(c(1, 0, 1))
#' fscore(X, X)
#' @export
fscore <- function(A, B, eps = 1e-4) {
  return(metrics(A, B, eps)[1])
}


# Compute the prial value between two matrices
# @param Ltrue true Laplacian matrix
# @param Lest estimated Laplacian matrix
# @param Lscm estimated Laplacian matrix via the generalized inverse of the
#        of the sample covariance matrix
prial <- function(Ltrue, Lest, Lscm) {
  return(100 * (1 - (norm(Lest - Ltrue, type = "F") /
                     norm(Lscm - Ltrue, type = "F"))^2))
}
