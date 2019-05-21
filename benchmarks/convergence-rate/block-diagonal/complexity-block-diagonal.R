library(igraph)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)

eps <- 1e-2
n_realizations <- 10
ratios <- c(30)
n <- c(20, 40, 80, 120, 160, 200, 400)
k <- 4
P <- diag(1, k)
maxiter <- 5e4
times <- matrix(0, n_realizations, length(n))

for (j in 1:length(n)) {
  mgraph <- sample_sbm(n[j], pref.matrix = P, block.sizes = c(rep(n[j] / k, k)))
  t <- as.integer(ratios * n[j])
  cat("\nRunning simulation for", t, "samples per node, t/n = ", ratios, "\n")
  for (r in 1:n_realizations) {
    print(r)
    E(mgraph)$weight <- runif(gsize(mgraph), min = 1e-1, max = 3)
    Ltrue <- as.matrix(laplacian_matrix(mgraph))
    # sample data from GP with covariance matrix set as
    # the pseudo inverse of the true Laplacian
    Y <- MASS::mvrnorm(t, mu = rep(0, n[j]), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    graph <- learn_k_component_graph(S, w0 = "qp", k = 4, beta = 1e2, fix_beta = TRUE,
                                     abstol = 0, maxiter = maxiter, verbose = TRUE)
    times[r, j] <- graph$elapsed_time[length(graph$elapsed_time)]
  }
}

saveRDS(times, file = "time-complexity.RDS")
