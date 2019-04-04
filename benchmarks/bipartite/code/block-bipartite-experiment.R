library(spectralGraphTopology)
library(igraph)
set.seed(234)

eps <- 5e-2
n1 <- 20
n2 <- 8
n <- 64
n_realizations <- 10
ratios <- c(5, 10, 30, 100, 250, 500, 1000)
rel_err_sgl <- matrix(0, n_realizations, length(ratios))
rel_err_naive <- matrix(0, n_realizations, length(ratios))
rel_err_qp <- matrix(0, n_realizations, length(ratios))
fscore_sgl <- matrix(0, n_realizations, length(ratios))
fscore_naive <- matrix(0, n_realizations, length(ratios))
fscore_qp <- matrix(0, n_realizations, length(ratios))
n_ratios <- c(1:length(ratios))

for (j in n_ratios) {
  p <- as.integer(ratios[j] * n)
  cat("\nRunning simulation for", p, "samples per node, p/n = ", ratios[j], "\n")
  for (r in 1:n_realizations) {
    b1 <- sample_bipartite(n1, n2, type="Gnp", p = .5, directed=FALSE)
    b2 <- sample_bipartite(n1-8, n2, type="Gnp", p = .6, directed=FALSE)
    b3 <- sample_bipartite(n1-12, n2, type="Gnp", p = .7, directed=FALSE)
    E(b1)$weight <- runif(gsize(b1), min = .1, max = 3)
    E(b2)$weight <- runif(gsize(b2), min = .1, max = 3)
    E(b3)$weight <- runif(gsize(b3), min = .1, max = 3)
    Lw1 <- as.matrix(laplacian_matrix(b1))
    Aw1 <- diag(diag(Lw1)) - Lw1
    Lw2 <- as.matrix(laplacian_matrix(b2))
    Aw2 <- diag(diag(Lw2)) - Lw2
    Lw3 <- as.matrix(laplacian_matrix(b3))
    Aw3 <- diag(diag(Lw3)) - Lw3
    Lw <- blockDiag(Lw1, Lw2, Lw3)
    Aw <- blockDiag(Aw1, Aw2, Aw3)
    # set number of samples
    Y <- MASS::mvrnorm(p, rep(0, n), Sigma = MASS::ginv(Lw))
    S <- cov(Y)
    Sinv <- MASS::ginv(S)
    w_naive <- spectralGraphTopology:::w_init(w0 = "naive", Sinv)
    w_qp <- spectralGraphTopology:::w_init(w0 = "qp", Sinv)
    graph <- learn_adjacency_and_laplacian(S, k = 3, z = 0, w0 = w_qp, beta = 1e3,
                                           fix_beta = TRUE, nu = 1e3, reltol = 1e-3, maxiter = 5e5)
    print(graph$convergence)
    Anaive <- A(w_naive)
    Aqp <- A(w_qp)
    rel_err_sgl[r, j] <- relativeError(Aw, graph$Adjacency)
    fscore_sgl[r, j] <- metrics(Aw, graph$Adjacency, eps)[1]
    rel_err_naive[r, j] <- relativeError(Aw, Anaive)
    fscore_naive[r, j] <- metrics(Aw, Anaive, eps)[1]
    rel_err_qp[r, j] <- relativeError(Aw, Aqp)
    fscore_qp[r, j] <- metrics(Aw, Aqp, eps)[1]
  }
  print(rel_err_sgl)
  print(rel_err_qp)
  print(fscore_sgl)
  print(fscore_qp)
}
saveRDS(rel_err_sgl, file = "block-rel-err-SGL.rds")
saveRDS(fscore_sgl, file = "block-fscore-SGL.rds")
saveRDS(rel_err_naive, file = "block-rel-err-naive.rds")
saveRDS(fscore_naive, file = "block-fscore-naive.rds")
saveRDS(rel_err_qp, file = "block-rel-err-QP.rds")
saveRDS(fscore_qp, file = "block-fscore-QP.rds")
