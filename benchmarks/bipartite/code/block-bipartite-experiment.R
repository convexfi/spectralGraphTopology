library(spectralGraphTopology)
library(igraph)
set.seed(234)

n1 <- 10
n2 <- 4
n <- 3*(n1+n2)-10
n_realizations <- 10
ratios <- c(5, 10, 30, 100, 250, 500, 1000)
rel_err_sgl <- array(0, length(ratios))
rel_err_naive <- array(0, length(ratios))
rel_err_qp <- array(0, length(ratios))
fscore_sgl <- array(0, length(ratios))
fscore_naive <- array(0, length(ratios))
fscore_qp <- array(0, length(ratios))
n_ratios <- c(1:length(ratios))

for (j in n_ratios) {
  p <- as.integer(ratios[j] * n)
  cat("\nRunning simulation for", p, "samples per node, p/n = ", ratios[j], "\n")
  for (r in 1:n_realizations) {
    b1 <- sample_bipartite(n1, n2, type="Gnp", p = .7, directed=FALSE)
    b2 <- sample_bipartite(n1-4, n2, type="Gnp", p = .8, directed=FALSE)
    b3 <- sample_bipartite(n1-6, n2, type="Gnp", p = .9, directed=FALSE)
    E(b1)$weight <- runif(gsize(b1), min = 1, max = 3)
    E(b2)$weight <- runif(gsize(b2), min = 1, max = 3)
    E(b3)$weight <- runif(gsize(b3), min = 1, max = 3)
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
    graph <- learn_adjacency_and_laplacian(S, k = 3, z = 0,
                                           w0 = w_qp, beta1 = 1e3, beta2 = 1e3, ftol = 1e-4, maxiter = 1e5)
    print(graph$convergence)
    Anaive <- A(w_naive)
    Aqp <- A(w_qp)
    rel_sgl = relativeError(Aw, graph$Aw)
    fs_sgl = Fscore(Aw, graph$Aw, 1e-1)
    rel_naive = relativeError(Aw, Anaive)
    fs_naive = Fscore(Aw, Anaive, 1e-1)
    rel_qp = relativeError(Aw, Aqp)
    fs_qp = Fscore(Aw, Aqp, 1e-1)
    rel_err_sgl[j] <- rel_err_sgl[j] + rel_sgl
    fscore_sgl[j] <- fscore_sgl[j] + fs_sgl
    rel_err_naive[j] <- rel_err_naive[j] + rel_naive
    fscore_naive[j] <- fscore_naive[j] + fs_naive
    rel_err_qp[j] <- rel_err_qp[j] + rel_qp
    fscore_qp[j] <- fscore_qp[j] + fs_qp
  }
  rel_err_sgl[j] <- rel_err_sgl[j] / n_realizations
  fscore_sgl[j] <- fscore_sgl[j] / n_realizations
  rel_err_naive[j] <- rel_err_naive[j] / n_realizations
  fscore_naive[j] <- fscore_naive[j] / n_realizations
  rel_err_qp[j] <- rel_err_qp[j] / n_realizations
  fscore_qp[j] <- fscore_qp[j] / n_realizations
  cat("\n** spectralGraphTopology results **\n")
  cat("Avg Relative error: ")
  cat(rel_err_sgl[j], "\n")
  cat("Avg Fscore: ")
  cat(fscore_sgl[j], "\n")
  cat("\n** Naive results **\n")
  cat("Avg Relative error: ")
  cat(rel_err_naive[j], "\n")
  cat("Avg Fscore: ")
  cat(fscore_naive[j], "\n")
  cat("\n** QP results **\n")
  cat("Avg Relative error: ")
  cat(rel_err_qp[j], "\n")
  cat("Avg Fscore: ")
  cat(fscore_qp[j], "\n")
}
saveRDS(rel_err_sgl, file = "block-rel-err-SGL.rds")
saveRDS(fscore_sgl, file = "block-fscore-SGL.rds")
saveRDS(rel_err_naive, file = "block-rel-err-naive.rds")
saveRDS(fscore_naive, file = "block-fscore-naive.rds")
saveRDS(rel_err_qp, file = "block-rel-err-QP.rds")
saveRDS(fscore_qp, file = "block-fscore-QP.rds")
