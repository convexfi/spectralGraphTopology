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


bipartite.obj_fun <- function(Aw, Lw, V, psi, K, J, nu) {
  return(bipartite.likelihood(Lw = Lw, K = K, J = J) +
         bipartite.prior(nu = nu, Aw = Aw, psi = psi, V = V))
}


bipartite.likelihood <- function(Lw, K, J) {
  return(sum(-log(eigval_sym(Lw + J)) + c(diag(K %*% Lw))))
}


bipartite.prior <- function(nu, Aw, psi, V) {
  return(.5 * nu * norm(Aw - V %*% diag(psi) %*% t(V), type="F")^2)
}


joint.obj_fun <- function(Lw, Aw, U, V, lambda, psi, beta, nu, K) {
  return(joint.likelihood(Lw = Lw, lambda = lambda, K = K) +
         joint.prior(beta = beta, nu = nu, Lw = Lw, Aw = Aw, U = U, V = V,
                     lambda = lambda, psi = psi))
}


joint.likelihood <- function(...) {
  return(laplacian.likelihood(...))
}


joint.prior <- function(beta, nu, Lw, Aw, U, V, lambda, psi) {
  return(laplacian.prior(beta = beta, Lw = Lw, lambda = lambda, U = U) +
         bipartite.prior(nu = nu, Aw = Aw, psi = psi, V = V))
}


dregular.obj_fun <- function(Lw, Aw, U, lambda, beta, eta, K, d) {
  return(dregular.likelihood(Lw, lambda, K) +
         dregular.prior(beta1, beta2, Lw, Aw, U, lambda, d))
}


dregular.likelihood <- function(...) {
  return(laplacian.likelihood(...))
}


dregular.prior <- function(beta, eta, Lw, Aw, U, lambda, d) {
  return(laplacian.prior(beta, Lw, lambda, U) +
         .5 * eta * norm(Aw + Lw - diag(d, ncol(Lw)), type="F")^2)
}


normalized_laplacian.lagragian <- function(lambda, Theta, K, M, U, beta) {
  ULmdUT <- crossprod(sqrt(lambda) * t(U))
  return(- sum(log(lambda)) + sum(diag(Theta %*% K))
         + sum(diag(t(M) %*% (Theta - crossprod(sqrt(lambda) * t(U)))))
         + .5 * beta * norm(Theta - ULmdUT, type="F")^2)
}
