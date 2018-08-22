#' Constructs a block diagonal matrix from a list of square matrices
#'
#' @param ... list of matrices or individual matrices
#' @return block diagonal matrix
blockDiag <- function(...) {
  return(blockDiagCpp(list(...)))
}


relativeError <- function(Xtrue, Xest) {
  return (100 * norm(Xtrue - Xest, type = "F") / max(1., norm(Xtrue, type = "F")))
}


prial <- function(Xtrue, Xest) {
  Xnaive <- MASS::ginv(cov(Xtrue))
  return (100 * (1 - (norm(Xest - Xtrue, type = "F") /
                      norm(Xnaive - Xtrue, type = "F"))^2))
}

# naive estimator
naive <- function(S) {
  return(MASS::ginv(S))
}

# qp estimator
qp <- function(S) {
  Sinv <- MASS::ginv(S)
  R <- vecLmat(ncol(Sinv))
  qp <- quadprog::solve.QP(t(R) %*% R, t(R) %*% vec(Sinv), diag(ncol(R)))
  return (L(qp$solution))
}
