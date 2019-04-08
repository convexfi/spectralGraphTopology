library(igraph)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)

set.seed(42)
N_realizations <- 5
ratios <- c(10, 100, 1000)
n_ratios <- c(1:length(ratios))
beta_set <- c(1, 10, 20, 50, 1e2, 1e3, 1e4, 1e5)
beta_size <- c(1:length(beta_set))
beta_rel_err <- matrix(0, length(beta_set), length(ratios))
beta_fscore <- matrix(0, length(beta_set), length(ratios))

n1 <- 10
n2 <- 4
n <- 3 * n1 + 3 * n2 - 10
b1 <- sample_bipartite(n1, n2, type="Gnp", p = .7, directed=FALSE)
b2 <- sample_bipartite(n1-4, n2, type="Gnp", p = .8, directed=FALSE)
b3 <- sample_bipartite(n1-6, n2, type="Gnp", p = .9, directed=FALSE)
for (j in n_ratios) {
  T <- as.integer(ratios[j] * n)
  cat("\nRunning simulation for", T, "samples per node, n / p = ", ratios[j], "\n")
  for (nrel in 1:N_realizations) {
    E(b1)$weight <- runif(gsize(b1), min = 0, max = 1)
    E(b2)$weight <- runif(gsize(b2), min = 0, max = 1)
    E(b3)$weight <- runif(gsize(b3), min = 0, max = 1)
    Lw1 <- as.matrix(laplacian_matrix(b1))
    Lw2 <- as.matrix(laplacian_matrix(b2))
    Lw3 <- as.matrix(laplacian_matrix(b3))
    Aw1 <- diag(diag(Lw1)) - Lw1
    Aw2 <- diag(diag(Lw2)) - Lw2
    Aw3 <- diag(diag(Lw3)) - Lw3
    Ltrue <- blockDiag(Lw1, Lw2, Lw3)
    Atrue <- blockDiag(Aw1, Aw2, Aw3)
    Y <- MASS::mvrnorm(T, mu = rep(0, n), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    for (i in beta_size) {
      print(beta_set[i])
      graph <- learn_adjacency_and_laplacian(S, w0 = "qp", k = 4, beta = beta_set[i], nu = beta_set[i],
                                             fix_beta = TRUE, maxiter = 100000)
      beta_rel_err[i, j] <- beta_rel_err[i, j] + relativeError(Ltrue, graph$Laplacian)
      beta_fscore[i, j] <- beta_fscore[i, j] + metrics(Ltrue, graph$Laplacian, 5e-2)[1]
    }
  }
  beta_rel_err[,j] <- beta_rel_err[,j] / N_realizations
  beta_fscore[,j] <- beta_fscore[,j] / N_realizations
}

saveRDS(beta_rel_err, file = "beta_gamma_rel_err.rds")
saveRDS(beta_fscore, file = "beta_gamma_fscore.rds")
