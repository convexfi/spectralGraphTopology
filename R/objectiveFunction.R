laplacian.objectiveFunction <- function(Lw, U, lambda, K, beta) {
  return(laplacian.likelihood(Lw, lambda, K) +
         laplacian.prior(beta, Lw, lambda, U))
}

laplacian.likelihood <- function(Lw, lambda, K) {
  return(sum(-log(lambda)) + sum(diag(K %*% Lw)))
}

laplacian.prior <- function(beta, Lw, lambda, U) {
  return(.5 * beta * norm(Lw - crossprod(sqrt(lambda) * t(U)), type="F")^2)
}

bipartite.obj_fun <- function(Aw, Lw, V, psi, K, J, beta) {
  return(bipartite.likelihood(Lw, K, J) +
         bipartite.prior(beta, Aw, psi, V))
}

bipartite.likelihood <- function(Lw, K, J) {
  return(sum(-log(eigenvalues(Lw + J)) + c(diag(K %*% Lw))))
}

bipartite.prior <- function(beta, Aw, psi, V) {
  return(.5 * beta * norm(Aw - V %*% diag(psi) %*% t(V), type="F")^2)
}

joint.obj_fun <- function(Lw, Aw, U, V, lambda, psi, beta1, beta2, K) {
  return(joint.likelihood(Lw, lambda, K) +
         joint.prior(beta1, beta2, Lw, Aw, U, V, lambda, psi))
}

joint.likelihood <- function(...) {
  return(laplacian.likelihood(...))
}

joint.prior <- function(beta1, beta2, Lw, Aw, U, V, lambda, psi) {
  return(laplacian.prior(beta1, Lw, lambda, U) + bipartite.prior(beta2, Aw, psi, V))
}

dregular.obj_fun <- function(Lw, Aw, U, lambda, beta1, beta2, K) {
  return(dregular.likelihood(Lw, lambda, K) +
         dregular.prior(beta1, beta2, Lw, Aw, U, lambda))
}

dregular.likelihood <- function(...) {
  return(laplacian.likelihood(...))
}

dregular.prior <- function(beta1, beta2, Lw, Aw, U, lambda) {
  return(laplacian.prior(beta1, Lw, lambda, U) +
         .5 * beta2 * norm(Aw - U %*% diag(d - lambda) %*% t(U), type="F")^2)
}
