# Negative log of the posterior density function
laplacian.objectiveFunction <- function(Lw, U, lambda, Kmat, beta) {
  return(laplacian.logLikelihood(Lw, lambda, Kmat) +
         laplacian.logPrior(beta, Lw, lambda, U))
}

# Negative logarithm of the likelihood function
laplacian.logLikelihood <- function(Lw, lambda, Kmat) {
  return(sum(-log(lambda)) + sum(diag(Kmat %*% Lw)))
}

# Negative logarithm of the prior probability distribution
laplacian.logPrior <- function(beta, Lw, lambda, U) {
  return(.5 * beta * norm(Lw - crossprod(sqrt(lambda) * t(U)), type="F")^2)
}

adjacency.objectiveFunction <- function(Aw, Lw, V, psi, Kmat, J, beta) {
  return(adjacency.logLikelihood(Lw, Kmat, J) +
         adjacency.logPrior(beta, Aw, psi, V))
}

adjacency.logLikelihood <- function(Lw, Kmat, J) {
  return(sum(-log(eigenvalues(Lw + J)) + c(diag(Kmat %*% Lw))))
}

adjacency.logPrior <- function(beta, Aw, psi, V) {
  return(.5 * beta * norm(Aw - V %*% diag(psi) %*% t(V), type="F")^2)
}

joint.objectiveFunction <- function(Lw, Aw, U, V, lambda, psi, beta1, beta2, Kmat) {
  return(joint.logLikelihood(Lw, lambda, Kmat) +
         joint.logPrior(beta1, beta2, Lw, Aw, U, V, lambda, psi))
}

joint.logLikelihood <- function(Lw, lambda, Kmat) {
  return(laplacian.logLikelihood(Lw, lambda, Kmat))
}

joint.logPrior <- function(beta1, beta2, Lw, Aw, U, V, lambda, psi) {
  return(.5 * beta1 * norm(Lw - crossprod(sqrt(lambda) * t(U))) +
         .5 * beta2 * norm(Aw - V %*% diag(psi) %*% t(V)))
}
