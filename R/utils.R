#' Constructs a block diagonal matrix from a list of square matrices
#'
#' @param ... list of matrices or individual matrices
#' @return block diagonal matrix
#' @export
blockDiag <- function(...) {
  return(blockDiagCpp(list(...)))
}

#' Compute the relative error between two matrices
#' @export
relativeError <- function(Xtrue, Xest) {
  return (norm(Xtrue - Xest, type = "F") / max(1., norm(Xtrue, type = "F")))
}


#' Compute the prial value between two matrices
#' @export
prial <- function(Xtrue, Xest) {
  Xnaive <- MASS::ginv(cov(Xtrue))
  return (100 * (1 - (norm(Xest - Xtrue, type = "F") /
                      norm(Xnaive - Xtrue, type = "F"))^2))
}
