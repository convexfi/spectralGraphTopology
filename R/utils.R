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


#' Computes the relative error between the true and estimated matrices
#'
#' @param West estimated matrix
#' @param Wtrue true matrix
#' @examples
#' library(spectralGraphTopology)
#' X <- L(c(1, 0, 1))
#' relative_error(X, X)
#' @export
relative_error <- function(West, Wtrue) {
  return(norm(West - Wtrue, type = "F") / norm(Wtrue, type = "F"))
}


#' Computes the fscore between two matrices
#'
#' @param Wtrue true matrix
#' @param West estimated matrix
#' @param eps real number such that edges whose values are smaller than eps are
#'            not considered in the computation of the fscore
#' @examples
#' library(spectralGraphTopology)
#' X <- L(c(1, 0, 1))
#' fscore(X, X)
#' @export
fscore <- function(Wtrue, West, eps = 1e-4) {
  return(metrics(Wtrue, West, eps)[1])
}

#' Computes the recall between two matrices
#'
#' @param Wtrue true matrix
#' @param West estimated matrix
#' @param eps real number such that edges whose values are smaller than eps are
#'            not considered in the computation of the fscore
#' @examples
#' library(spectralGraphTopology)
#' X <- L(c(1, 0, 1))
#' recall(X, X)
#' @export
recall <- function(Wtrue, West, eps = 1e-4) {
  return(metrics(Wtrue, West, eps)[2])
}


#' Computes the specificity between two matrices
#'
#' @param Wtrue true matrix
#' @param West estimated matrix
#' @param eps real number such that edges whose values are smaller than eps are
#'            not considered in the computation of the fscore
#' @examples
#' library(spectralGraphTopology)
#' X <- L(c(1, 0, 1))
#' specificity(X, X)
#' @export
specificity <- function(Wtrue, West, eps = 1e-4) {
  return(metrics(Wtrue, West, eps)[3])
}



#' Computes the false discovery rate between two matrices
#'
#' @param Wtrue true matrix
#' @param West estimated matrix
#' @param eps real number such that edges whose values are smaller than eps are
#'            not considered in the computation of the fscore
#' @examples
#' library(spectralGraphTopology)
#' X <- L(c(1, 0, 1))
#' fdr(X, X)
#' @export
fdr <- function(Wtrue, West, eps = 1e-4) {
  return(metrics(Wtrue, West, eps)[6])
}


#' Computes the negative predictive value between two matrices
#'
#' @param Wtrue true matrix
#' @param West estimated matrix
#' @param eps real number such that edges whose values are smaller than eps are
#'            not considered in the computation of the fscore
#' @examples
#' library(spectralGraphTopology)
#' X <- L(c(1, 0, 1))
#' npv(X, X)
#' @export
npv <- function(Wtrue, West, eps = 1e-4) {
  return(metrics(Wtrue, West, eps)[5])
}


#' Computes the accuracy between two matrices
#'
#' @param Wtrue true matrix
#' @param West estimated matrix
#' @param eps real number such that edges whose values are smaller than eps are
#'            not considered in the computation of the fscore
#' @examples
#' library(spectralGraphTopology)
#' X <- L(c(1, 0, 1))
#' accuracy(X, X)
#' @export
accuracy <- function(Wtrue, West, eps = 1e-4) {
  return(metrics(Wtrue, West, eps)[3])
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
