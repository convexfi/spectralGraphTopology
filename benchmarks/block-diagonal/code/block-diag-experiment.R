library(igraph)
library(R.matlab)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)

set.seed(42)

eps <- 5e-2
N_realizations <- 5
ratios <- c(5e3, 1e3, 5e2, 1e2, 1e1, 1.5, .5)
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
P <- diag(.3, k)
# K-component graph

for (j in n_ratios) {
  T <- as.integer(ratios[j] * n)
  cat("\nRunning simulation for", T, "samples per node, p/n = ", ratios[j], "\n")
  for (i in 1:N_realizations) {
    mgraph <- sample_sbm(n, pref.matrix = P, block.sizes = c(rep(n / k, k)))
    E(mgraph)$weight <- runif(gsize(mgraph), min = 1e-1, max = 3)
    Ltrue <- as.matrix(laplacian_matrix(mgraph))
    Y <- MASS::mvrnorm(T, mu = rep(0, n), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    Lnaive <- MASS::ginv(S)
    w_qp <- spectralGraphTopology:::w_init("qp", Lnaive)
    Lqp <- L(w_qp)
    graph <- learn_laplacian_matrix(S, w0 = w_qp, k = 4, ub = 32, beta = 1e2, fix_beta = TRUE,
                                    edge_tol = 1e-3, abstol = 1e-4, maxiter = 5e5)
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
