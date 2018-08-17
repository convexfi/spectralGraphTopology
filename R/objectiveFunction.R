objectiveFunction <- function(Theta, U, lambda, Km, beta, N, K) {
  return(sum(-log(lambda)) + sum(diag(Km %*% Theta)) +
         .5 * beta * norm(Theta - crossprod(sqrt(lambda) * t(U)), type="F")^2)
}
