w_init <- function(w0, Sinv) {
  if (is.character(w0)) {
    if (w0 == "qp") {
      R <- vecLmat(ncol(Sinv))
      qp <- quadprog::solve.QP(crossprod(R), t(R) %*% vec(Sinv), diag(ncol(R)))
      w0 <- qp$solution
    } else if (w0 == "naive") {
      w0 <- Linv(Sinv)
    }
  }
  return(w0)
}


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
  c1 <- Lstar(beta * ULmdUT - K)
  c2 <- nu * Astar(VPsiVT)
  Mw <- Lstar(Lw)
  Pw <- 2 * w
  grad_f1 <- beta * Mw - c1
  M_grad_f1 <- Lstar(L(grad_f1))
  grad_f2 <- nu * Pw - c2
  P_grad_f2 <- 2 * grad_f2
  grad_f <- grad_f1 + grad_f2
  t <- sum((beta * Mw + nu * Pw - (c1 + c2)) * grad_f) / sum(grad_f * (beta * M_grad_f1 + nu * P_grad_f2))
  w_update <- w - t * (grad_f1 + grad_f2)
  w_update[w_update < 0] <- 0
  return(w_update)
}


bipartite.w_update <- function(w, Aw, V, nu, psi, K, J, Lips) {
  grad_h <- 2 * w - Astar(V %*% diag(psi) %*% t(V)) #+ Lstar(K) / beta
  w_update <- w - (Lstar(inv_sympd(L(w) + J) + K) + nu * grad_h) / (2 * nu + Lips)
  w_update[w_update < 0] <- 0
  return(w_update)
}


laplacian.U_update <- function(Lw, k) {
  return(eigvec_sym(Lw)[, (k+1):ncol(Lw)])
}


bipartite.V_update <- function(Aw, z) {
  n <- ncol(Aw)
  V <- eigvec_sym(Aw)
  return(cbind(V[, 1:(.5*(n - z))], V[, (1 + .5*(n + z)):n]))
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
  x <- stats::isoreg(rev(c_tilde))$yf
  x <- c(-rev(x), x)
  x[x < lb] = lb
  x[x > ub] = ub
  return(x)
}


joint.lambda_update <- function(...) {
  return(laplacian.lambda_update(...))
}


joint.psi_update <- function(...) {
  return(bipartite.psi_update(...))
}
