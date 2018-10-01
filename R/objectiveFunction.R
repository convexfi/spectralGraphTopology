# Negative log of the posterior density function
objectiveFunction <- function(Lw, U, lambda, Kmat, beta) {
  return(logLikelihood(Lw, lambda, Kmat) + logPrior(beta, Lw, lambda, U))
}

# Negative logarithm of the likelihood function
logLikelihood <- function(Lw, lambda, Kmat) {
  return(sum(-log(lambda)) + sum(diag(Kmat %*% Lw)))
}

# Negative logarithm of the prior probability distribution
logPrior <- function(beta, Lw, lambda, U) {
  return(.5 * beta * norm(Lw - crossprod(sqrt(lambda) * t(U)), type="F")^2)
}
