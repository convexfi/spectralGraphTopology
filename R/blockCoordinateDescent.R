#' Updates the value of w, the weight vector of the graph
#'
#' @param w weight vector of the graph
#' @param U matrix whose columns represent the eigen vectors of the Laplacian
#' matrix in increasing order
#' @param β scalar that controls the strength of the regularization term
#' @param λ vector whose entries are the eigenvalues of the Laplacian in
#' increasing order
#' @param N dimension of each data sample
#' @param Km matrix (see section 1.1 for its definition)
#'
#' @return the updated value of w
w_update <- function(w, U, beta, lambda, N, Km) {
  #grad_f <- Lstar(L(w) - U %*% diag(lambda) %*% t(U) + Km / beta)
  grad_f <- Lstar(L(w) - crossprod(sqrt(lambda) * t(U)) + Km / beta)
  w_update <- w - .5 * grad_f / N
  return(pmax(0, w_update))
}


#' Updates the value of U, the matrix whose columns represent the eigen vectors
#' of the Laplacian in increasing order
#'
#' @param w weight vector of the graph
#' @param N dimension of each data sample
#'
#' @return the updated value of U
U_update <- function(w, N) {
  return(eigen(L(w))$vectors[, N:1])
}


#' Updates the value of lambda, the vector whose entries correspond to the
#' eigenvalues of the Laplacian matrix
#'
#' @param lb lower bound on the eigenvalues of the Laplacian matrix
#' @param ub upper bound on the eigenvalues of the Laplacian matrix
#' @param β: scalar that controls the strength of the regularization term
#' @param U matrix whose columns represent the eigen vectors of the Laplacian
#' matrix in increasing order
#' @param w weight vector of the graph
#' @param N dimension of each data sample
#' @param K number of components of the graph
#'
#' @return the updated value of λ
lambda_update <- function(lb, ub, beta, U, w, N, K) {
  d <- diag(t(U) %*% L(w) %*% U)
  lambda <- .5 * (d + sqrt(d^2 + 4 / beta)) # unconstrained solution as initial point
  lambda[1:K] <- 0.
  condition <- all(lambda[N] <= ub, lambda[K+1] >= lb,
                   lambda[(K+2):N] >= lambda[(K+1):(N-1)])

  if (condition) {
    return (lambda)
  } else {
    bigger_ub <- lambda[(K+1):N] > ub
    smaller_lb <- lambda[(K+1):N] < lb
    lambda[(K+1):N][bigger_ub] <- ub
    lambda[(K+1):N][smaller_lb] <- lb
  }

  condition <- all(lambda[N] <= ub, lambda[K+1] >= lb,
                   lambda[(K+2):N] >= lambda[(K+1):(N-1)])

  if (condition) {
    return (lambda)
  } else {
    stop("eigenvalues are not in increasing order")
  }
}
