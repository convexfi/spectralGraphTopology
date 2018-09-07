#' Given a reference Laplacian matrix, this function computes the relative error,
#' the negative loglikelihood, and the objective function per iteration for the
#' spectralGraphTopology algorithm
#'
#' @param Lw reference Laplacian matrix from which random samples will be drawn
#' @param T number of samples to be drawn which will be used to estimate L
#' @param K number of components of the graph
#' @param ... extra arguments to be passed to learnGraphTopology
#' @return rel_err sequence of relative error measurements for each iteration
#' @return obj_fun objective function value for each iteration
#' @return loglike negative loglikelihood value for each iteration
#' @export
performance_per_iteration <- function(Lw, T, K, ...) {
  N <- ncol(Lw)
  Y <- MASS::mvrnorm(T, rep(0, N), MASS::ginv(Lw))
  lgt <- learnGraphTopology(Y, K, ...)
  relerr_seq <- c()
  ii <- c(1:length(lgt$w_seq))
  for (i in ii) {
    relerr <- relativeError(Lw, L(lgt$w_seq[[i]]))
    relerr_seq <- c(relerr_seq, relerr)
  }
  return(list(rel_err = relerr_seq, obj_fun = lgt$obj_fun,
              loglike = lgt$loglike))
}
