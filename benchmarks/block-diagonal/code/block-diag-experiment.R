library(igraph)
library(R.matlab)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)

set.seed(42)

eps <- 5e-2
N_realizations <- 5
ratios <- c(10, 50, 100, 5e2, 1e3)
n_ratios <- c(1:length(ratios))

rel_err_sgl <- matrix(0, N_realizations, length(ratios))
rel_err_naive <- matrix(0, N_realizations, length(ratios))
rel_err_qp <- matrix(0, N_realizations, length(ratios))
recall_sgl <- matrix(0, N_realizations, length(ratios))
recall_naive <- matrix(0, N_realizations, length(ratios))
recall_qp <- matrix(0, N_realizations, length(ratios))
specificity_sgl <- matrix(0, N_realizations, length(ratios))
specificity_naive <- matrix(0, N_realizations, length(ratios))
specificity_qp <- matrix(0, N_realizations, length(ratios))
accuracy_sgl <- matrix(0, N_realizations, length(ratios))
accuracy_naive <- matrix(0, N_realizations, length(ratios))
accuracy_qp <- matrix(0, N_realizations, length(ratios))
fscore_sgl <- matrix(0, N_realizations, length(ratios))
fscore_naive <- matrix(0, N_realizations, length(ratios))
fscore_qp <- matrix(0, N_realizations, length(ratios))

n <- 64
k <- 4
P <- diag(.5, k)
# K-component graph
mgraph <- sample_sbm(n, pref.matrix = P, block.sizes = c(rep(n / k, k)))

for (j in n_ratios) {
  T <- as.integer(ratios[j] * n)
  cat("\nRunning simulation for", T, "samples per node, n/p = ", ratios[j], "\n")
  for (i in 1:N_realizations) {
    #graph <- learn_laplacian_matrix(S, w0 = w_qp, k = 4, ub = 32, beta = 20, fix_beta = TRUE,
    #                                abstol = 0, edge_tol = 0, maxiter = 5e5)

    E(mgraph)$weight <- runif(gsize(mgraph), min = 0, max = 1)
    Ltrue <- as.matrix(laplacian_matrix(mgraph))
    Y <- MASS::mvrnorm(T, mu = rep(0, n), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    graph <- learn_k_component_graph(S, w0 = "qp", k = k, beta = 4, fix_beta = TRUE,
                                    alpha = 1e-2, maxiter = 100000, abstol = 0)
    Sinv <- MASS::ginv(S)
    w_qp <- spectralGraphTopology:::w_init("qp", Sinv)
    w_naive <- spectralGraphTopology:::w_init("naive", Sinv)
    Lnaive <- L(w_naive)
    Lqp <- L(w_qp)
    metrics_sgl <- metrics(Ltrue, graph$Laplacian, eps)
    metrics_qp <- metrics(Ltrue, Lqp, eps)
    metrics_naive <- metrics(Ltrue, Lnaive, eps)
    rel_err_sgl[i, j] <- relativeError(Ltrue, graph$Laplacian)
    rel_err_qp[i, j] <- relativeError(Ltrue, Lqp)
    rel_err_naive[i, j] <- relativeError(Ltrue, Lnaive)
    fscore_sgl[i, j] <- metrics_sgl[1]
    fscore_qp[i, j] <- metrics_qp[1]
    fscore_naive[i, j] <- metrics_naive[1]
    recall_sgl[i, j] <- metrics_sgl[2]
    recall_qp[i, j] <- metrics_qp[2]
    recall_naive[i, j] <- metrics_naive[2]
    specificity_sgl[i, j] <- metrics_sgl[3]
    specificity_qp[i, j] <- metrics_qp[3]
    specificity_naive[i, j] <- metrics_naive[3]
    accuracy_sgl[i, j] <- metrics_sgl[4]
    accuracy_qp[i, j] <- metrics_qp[4]
    accuracy_naive[i, j] <- metrics_naive[4]
    print(rel_err_sgl)
    print(rel_err_qp)
    print(fscore_sgl)
    print(fscore_qp)
  }
}

saveRDS(rel_err_sgl, file = "rel-err-SGL.rds")
saveRDS(rel_err_naive, file = "rel-err-naive.rds")
saveRDS(rel_err_qp, file = "rel-err-QP.rds")
saveRDS(fscore_sgl, file = "fscore-SGL.rds")
saveRDS(fscore_naive, file = "fscore-naive.rds")
saveRDS(fscore_qp, file = "fscore-QP.rds")
saveRDS(recall_sgl, file = "recall-SGL.rds")
saveRDS(recall_naive, file = "recall-naive.rds")
saveRDS(recall_qp, file = "recall-QP.rds")
saveRDS(specificity_sgl, file = "specificity-SGL.rds")
saveRDS(specificity_naive, file = "specificity-naive.rds")
saveRDS(specificity_qp, file = "specificity-QP.rds")
saveRDS(accuracy_sgl, file = "accuracy-SGL.rds")
saveRDS(accuracy_naive, file = "accuracy-naive.rds")
saveRDS(accuracy_qp, file = "accuracy-QP.rds")
