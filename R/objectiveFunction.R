# Negative log of the posterior density function
laplacian.objectiveFunction <- function(Lw, U, lambda, Kmat, beta) {
  return(laplacian.loglikelihood(Lw, lambda, Kmat) +
         laplacian.logprior(beta, Lw, lambda, U))
}

# Negative loglikelihood function
laplacian.loglikelihood <- function(Lw, lambda, Kmat) {
  return(sum(-log(lambda)) + sum(diag(Kmat %*% Lw)))
}

# Negative log of the prior density function
laplacian.logprior <- function(beta, Lw, lambda, U) {
  return(.5 * beta * norm(Lw - crossprod(sqrt(lambda) * t(U)), type="F")^2)
}

bipartite.obj_fun <- function(Aw, Lw, V, psi, Kmat, J, beta) {
  return(bipartite.loglikelihood(Lw, Kmat, J) +
         bipartite.logprior(beta, Aw, psi, V))
}

bipartite.loglikelihood <- function(Lw, Kmat, J) {
  return(sum(-log(eigenvalues(Lw + J)) + c(diag(Kmat %*% Lw))))
}

bipartite.logprior <- function(beta, Aw, psi, V) {
  return(.5 * beta * norm(Aw - V %*% diag(psi) %*% t(V), type="F")^2)
}

joint.obj_fun <- function(Lw, Aw, U, V, lambda, psi, beta1, beta2, Kmat) {
  return(joint.loglikelihood(Lw, lambda, Kmat) +
         joint.logprior(beta1, beta2, Lw, Aw, U, V, lambda, psi))
}

joint.loglikelihood <- function(Lw, lambda, Kmat) {
  return(laplacian.loglikelihood(Lw, lambda, Kmat))
}

joint.logprior <- function(beta1, beta2, Lw, Aw, U, V, lambda, psi) {
  return(laplacian.logprior(beta1, Lw, lambda, U) + bipartite.logprior(beta2, Aw, psi, V))
}
