laplacian.objectiveFunction <- function(Lw, U, lambda, K, beta) {
  return(laplacian.loglikelihood(Lw, lambda, K) +
         laplacian.logprior(beta, Lw, lambda, U))
}

laplacian.loglikelihood <- function(Lw, lambda, K) {
  return(sum(-log(lambda)) + sum(diag(K %*% Lw)))
}

laplacian.logprior <- function(beta, Lw, lambda, U) {
  return(.5 * beta * norm(Lw - crossprod(sqrt(lambda) * t(U)), type="F")^2)
}

bipartite.obj_fun <- function(Aw, Lw, V, psi, K, J, beta) {
  return(bipartite.loglikelihood(Lw, K, J) +
         bipartite.logprior(beta, Aw, psi, V))
}

bipartite.loglikelihood <- function(Lw, K, J) {
  return(sum(-log(eigenvalues(Lw + J)) + c(diag(K %*% Lw))))
}

bipartite.logprior <- function(beta, Aw, psi, V) {
  return(.5 * beta * norm(Aw - V %*% diag(psi) %*% t(V), type="F")^2)
}

joint.obj_fun <- function(Lw, Aw, U, V, lambda, psi, beta1, beta2, K) {
  return(joint.loglikelihood(Lw, lambda, K) +
         joint.logprior(beta1, beta2, Lw, Aw, U, V, lambda, psi))
}

joint.loglikelihood <- function(...) {
  return(laplacian.loglikelihood(...))
}

joint.logprior <- function(beta1, beta2, Lw, Aw, U, V, lambda, psi) {
  return(laplacian.logprior(beta1, Lw, lambda, U) + bipartite.logprior(beta2, Aw, psi, V))
}

dregular.obj_fun <- function(Lw, Aw, U, lambda, beta1, beta2, K) {
  return(dregular.loglikelihood(Lw, lambda, K) +
         dregular.logprior(beta1, beta2, Lw, Aw, U, lambda))
}

dregular.loglikelihood <- function(...) {
  return(laplacian.loglikelihood(...))
}

dregular.logprior <- function(beta1, beta2, Lw, Aw, U, lambda) {
  return(laplacian.logprior(beta1, Lw, lambda, U) +
         .5 * beta2 * norm(Aw - U %*% diag(d - lambda) %*% t(U), type="F")^2)
}
