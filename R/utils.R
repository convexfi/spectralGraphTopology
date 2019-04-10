#' Constructs a block diagonal matrix from a list of square matrices
#'
#' @param ... list of matrices or individual matrices
#' @return block diagonal matrix
#' @export
block_diag <- function(...) {
  return(blockDiagCpp(list(...)))
}


#' Compute the relative error between two matrices
#' @param Ltrue true Laplacian matrix
#' @param Lest estimated Laplacian matrix
#' @export
relative_error <- function(Ltrue, Lest) {
  return(norm(Ltrue - Lest, type = "F") / norm(Ltrue, type = "F"))
}


#' Compute the prial value between two matrices
#' @param Ltrue true Laplacian matrix
#' @param Lest estimated Laplacian matrix
#' @param Lscm estimated Laplacian matrix via the generalized inverse of the
#'        of the sample covariance matrix
#' @export
prial <- function(Ltrue, Lest, Lscm) {
  return(100 * (1 - (norm(Lest - Ltrue, type = "F") /
                     norm(Lnaive - Ltrue, type = "F"))^2))
}
