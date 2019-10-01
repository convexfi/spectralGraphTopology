library(CVXR)

Theta <- Variable(p, p)
J <- matrix(1/p, p, p)
objective <- sum(diag(Theta %*% (K + J))) - log_det(Theta + J)
problem <- Minimize()
