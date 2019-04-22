library(igraph)
library(R.matlab)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)

set.seed(42)

eps <- 5e-2
N_realizations <- 5
ratios <- c(2, 5, 10, 50, 100, 5e2, 1e3)
n_ratios <- c(1:length(ratios))

accuracy <- matrix(0, N_realizations, length(ratios))
fscore <- matrix(0, N_realizations, length(ratios))

n <- 64
k <- 4
P <- diag(.5, k)
# K-component graph
mgraph <- sample_sbm(n, pref.matrix = P, block.sizes = c(rep(n / k, k)))

for (j in n_ratios) {
  T <- as.integer(ratios[j] * n)
  cat("\nRunning simulation for", T, "samples per node, n/p = ", ratios[j], "\n")
  for (i in 1:N_realizations) {
    E(mgraph)$weight <- runif(gsize(mgraph), min = 0, max = 1)
    Ltrue <- as.matrix(laplacian_matrix(mgraph))
    Y <- MASS::mvrnorm(T, mu = rep(0, n), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    graph <- learn_normalized_laplacian(S, k = 4, maxiter = 100000)
    metrics <- metrics(Ltrue, graph$Laplacian, eps)
    fscore[i, j] <- metrics[1]
    accuracy[i, j] <- metrics[4]
    print(fscore)
    print(accuracy)
  }
}

saveRDS(fscore, file = "fscore-NSGL.rds")
saveRDS(accuracy, file = "accuracy-NSGL.rds")
