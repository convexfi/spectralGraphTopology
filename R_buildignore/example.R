library(spectralGraphTopology)

# Number of columns of Lw
n <- 4
# Number of random samples to generate
k <- 100000
# True w vector
w <- runif(6)
# Number of components
K <- 1
# True Lw matrix
Lw <- L(w)
# Generate random samples
Y <- t(MASS::mvrnorm(k, rep(0, n), MASS::ginv(Lw)))
# Learn Lw
Lw_est <- learnGraphTopology(Y, K, lb=1e-4, ub=1e2, beta=1., rho=.5, maxiter=5000, Lw_tol = 1e-12)
print(Lw)
print(Lw_est)
print(norm(Lw - Lw_est, type="F") / max(1., norm(Lw, type="F")))

# Results should be something similiar to
#     [,1] [,2] [,3] [,4]
# [1,]    6   -1   -2   -3
# [2,]   -1   10   -4   -5
# [3,]   -2   -4   12   -6
# [4,]   -3   -5   -6   14
#     [,1]       [,2]      [,3]      [,4]
# [1,]  5.9054503 -0.9446548 -1.947063 -3.013733
# [2,] -0.9446548 10.0907062 -4.068343 -5.077708
# [3,] -1.9470629 -4.0683432 12.017785 -6.002379
# [4,] -3.0137326 -5.0777083 -6.002379 14.093820
# [1] 0.009535078




#w1 <- runif(6)
#w2 <- runif(6)
#Theta1 <- CppLOp(w1, n)
#Theta2 <- CppLOp(w2, n)
#Theta <- rbind(cbind(Theta1, matrix(0, n, n)),
#               cbind(matrix(0, n, n), Theta2))
#eigen(Theta)$values
