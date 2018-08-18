#' Updates the value of w, the weight vector of the graph
#'
#' @param w weight vector of the graph
#' @param U matrix whose columns represent the eigen vectors of the Laplacian
#' matrix in increasing order
#' @param β scalar that controls the strength of the regularization term
#' @param λ vector whose entries are the eigenvalues of the Laplacian in
#' increasing order
#' @param N dimension of each data sample
#' @param Kmat matrix (see section 1.1 for its definition)
#'
#' @return the updated value of w
w_update <- function(w, U, beta, lambda, N, Kmat) {
  #grad_f <- Lstar(L(w) - U %*% diag(lambda) %*% t(U) + Kmat / beta)
  grad_f <- Lstar(L(w) - crossprod(sqrt(lambda) * t(U)) + Kmat / beta)
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
U_update <- function(w, N, K) {
  return(eigen(L(w))$vectors[, (N-K):1])
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
  q <- N - K
  d <- diag(t(U) %*% L(w) %*% U)
  lambda <- .5 * (d + sqrt(d^2 + 4 / beta)) # unconstrained solution as initial point
  condition <- all(lambda[q] <= ub, lambda[1] >= lb,
                   lambda[2:q] >= lambda[1:(q-1)])

  if (condition) {
    return (lambda)
  } else {
    greater_ub <- lambda > ub
    lesser_lb <- lambda < lb
    lambda[greater_ub] <- ub
    lambda[lesser_lb] <- lb
  }

  condition <- all(lambda[q] <= ub, lambda[1] >= lb,
                   lambda[2:q] >= lambda[1:(q-1)])

  if (condition) {
    return (lambda)
  } else {
    stop("eigenvalues are not in increasing order")
  }
}
