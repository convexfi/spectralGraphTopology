library(spectralGraphTopology)

# Number of samples
T <- 100000
# True w vector
w <- sample(1:10, 10)
# Number of components
K <- 1
# True Theta matrix
Theta <- L(w)
# Generate random samples
N <- as.integer(.5 * (1 + sqrt(1 + 8 * length(w))))
Y <- t(MASS::mvrnorm(T, rep(0, N), MASS::ginv(Theta)))
# Learn Theta
Theta_est <- learnGraphTopology(Y, K, lb=1e-4, ub=1e2, beta=.1, maxiter=5000)
print(Theta)
print(Theta_est)
print(norm(Theta - Theta_est, type="F") / max(1., norm(Theta, type="F")))

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
