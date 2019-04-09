w_init <- function(w0, Sinv) {
  if (is.character(w0)) {
    if (w0 == "qp") {
      R <- Matrix::Matrix(vecLmat(ncol(Sinv)), sparse = TRUE)
      qp <- osqp::solve_osqp(P = Matrix::crossprod(R), q = -Matrix::t(R) %*% vec(Sinv),
                             A = Matrix::Matrix(diag(ncol(R)), sparse = TRUE),
                             l = rep_len(0, ncol(R)), u = rep_len(Inf, ncol(R)),
                             osqp::osqpSettings(verbose = FALSE))
      w0 <- qp$x
    } else if (w0 == "naive") {
      w0 <- Linv(Sinv)
      w0[w0 < 0] <- 0
    }
  }
  return(w0)
}


#laplacian.w_update <- function(w, Lw, U, beta, lambda, K) {
#  grad_f <- Lstar(Lw - crossprod(sqrt(lambda) * t(U)) + K / beta)
#  w_update <- w - .5 * grad_f / nrow(Lw)
#  w_update[w_update < 0] <- 0
#  return(w_update)
#}


laplacian.w_update <- function(w, Lw, U, beta, lambda, K) {
  c <- Lstar(crossprod(sqrt(lambda) * t(U)) - K / beta)
  grad_f <- Lstar(Lw) - c
  M_grad_f <- Lstar(L(grad_f))
  wT_M_grad_f <- sum(w * M_grad_f)
  dwT_M_dw <- sum(grad_f * M_grad_f)
  # exact line search
  t <- (wT_M_grad_f - sum(c * grad_f)) / dwT_M_dw
  w_update <- w - t * grad_f
  w_update[w_update < 0] <- 0
  return(w_update)
}


joint.w_update <- function(w, Lw, Aw, U, V, lambda, psi, beta, nu, K) {
  ULmdUT <- crossprod(sqrt(lambda) * t(U))
  VPsiVT <- V %*% diag(psi) %*% t(V)
  grad_f1 <- Lstar(beta * (Lw - ULmdUT) + K)
  grad_f2 <- nu * Astar(Aw - VPsiVT)
  w_update <- w - .5 * (grad_f1 + grad_f2) / (nrow(Lw) * beta + nu)
  w_update[w_update < 0] <- 0
  return(w_update)
}


bipartite.w_update <- function(w, Aw, V, nu, psi, K, J, Lips) {
  grad_h <- 2 * w - Astar(V %*% diag(psi) %*% t(V)) #+ Lstar(K) / beta
  w_update <- w - (Lstar(inv_sympd(L(w) + J) + K) + nu * grad_h) / (2 * nu + Lips)
  w_update[w_update < 0] <- 0
  return(w_update)
}


dregular.w_update <- function(w, Lw, Aw, U, beta, eta, lambda, d, K) {
  n <- ncol(Aw)
  ULmdUT <- crossprod(sqrt(lambda) * t(U))
  grad_f <- Lstar(beta * (Lw - ULmdUT) + K + eta * diag(diag(Lw) - rep(d, n)))
  w_update <- w - .5 * grad_f / ((beta + eta) * n + eta)
  w_update[w_update < 0] <- 0
  return(w_update)
}


laplacian.U_update <- function(Lw, k) {
  return(eigenvectors(Lw)[, (k+1):ncol(Lw)])
}


bipartite.V_update <- function(Aw, z) {
  n <- ncol(Aw)
  V <- eigenvectors(Aw)
  return(cbind(V[, 1:(.5*(n - z))], V[, (1 + .5*(n + z)):n]))
}


dregular.U_update <- function(...) {
  return(laplacian.U_update(...))
}


joint.U_update <- function(...) {
  return(laplacian.U_update(...))
}


joint.V_update <- function(...) {
  return(bipartite.V_update(...))
}


laplacian.lambda_update <- function(lb, ub, beta, U, Lw, k) {
  q <- ncol(Lw) - k
  d <- diag(t(U) %*% Lw %*% U)
  # unconstrained solution as initial point
  lambda <- .5 * (d + sqrt(d^2 + 4 / beta))
  eps <- 1e-9
  condition <- c((lambda[q] - ub) <= eps,
                 (lambda[1] - lb) >= -eps,
                 (lambda[2:q] - lambda[1:(q-1)]) >= -eps)
  if (all(condition)) {
    return (lambda)
  } else {
    greater_ub <- lambda > ub
    lesser_lb <- lambda < lb
    lambda[greater_ub] <- ub
    lambda[lesser_lb] <- lb
  }
  condition <- c((lambda[q] - ub) <= eps,
                 (lambda[1] - lb) >= -eps,
                 (lambda[2:q] - lambda[1:(q-1)]) >= -eps)
  if (all(condition)) {
    return (lambda)
  } else {
    print(lambda)
    stop("eigenvalues are not in increasing order,
          consider increasing the value of beta")
  }
}


bipartite.psi_update <- function(V, Aw, lb = -Inf, ub = Inf) {
  c <- diag(t(V) %*% Aw %*% V)
  n <- length(c)
  c_tilde <- .5 * (rev(c[(n/2 + 1):n]) - c[1:(n/2)])
  x <- isoreg(rev(c_tilde))$yf
  x <- c(-rev(x), x)
  x[x < lb] = lb
  x[x > ub] = ub
  return(x)
}


dregular.lambda_update <- function(...) {
  return(laplacian.lambda_update(...))
}


joint.lambda_update <- function(...) {
  return(laplacian.lambda_update(...))
}


joint.psi_update <- function(...) {
  return(bipartite.psi_update(...))
}
