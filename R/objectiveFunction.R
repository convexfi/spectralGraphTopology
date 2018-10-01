# Negative log of the posterior density function
objectiveFunction <- function(Theta, U, lambda, Kmat, beta) {
  return(logLikelihood(Theta, lambda, Kmat) + logPrior(beta, Theta, lambda, U))
}

# Negative logarithm of the likelihood function
logLikelihood <- function(Theta, lambda, Kmat) {
  return(sum(-log(lambda)) + sum(diag(Kmat %*% Theta)))
}

# Negative logarithm of the prior probability distribution
logPrior <- function(beta, Theta, lambda, U) {
  return(.5 * beta * norm(Theta - crossprod(sqrt(lambda) * t(U)), type="F")^2)
}
