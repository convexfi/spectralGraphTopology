library(igraph)
library(R.matlab)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)

set.seed(42)

eps <- 5e-2
N_realizations <- 5
ratios <- rev(c(.5, 2, 10, 50, 100, 5e2, 1e3))
n_ratios <- c(1:length(ratios))

rel_err_sgl <- matrix(0, N_realizations, length(ratios))
rel_err_cgl <- matrix(0, N_realizations, length(ratios))
rel_err_naive <- matrix(0, N_realizations, length(ratios))
rel_err_qp <- matrix(0, N_realizations, length(ratios))
recall_sgl <- matrix(0, N_realizations, length(ratios))
recall_cgl <- matrix(0, N_realizations, length(ratios))
recall_naive <- matrix(0, N_realizations, length(ratios))
recall_qp <- matrix(0, N_realizations, length(ratios))
specificity_sgl <- matrix(0, N_realizations, length(ratios))
specificity_cgl <- matrix(0, N_realizations, length(ratios))
specificity_naive <- matrix(0, N_realizations, length(ratios))
specificity_qp <- matrix(0, N_realizations, length(ratios))
accuracy_sgl <- matrix(0, N_realizations, length(ratios))
accuracy_cgl <- matrix(0, N_realizations, length(ratios))
accuracy_naive <- matrix(0, N_realizations, length(ratios))
accuracy_qp <- matrix(0, N_realizations, length(ratios))
fscore_sgl <- matrix(0, N_realizations, length(ratios))
fscore_cgl <- matrix(0, N_realizations, length(ratios))
fscore_naive <- matrix(0, N_realizations, length(ratios))
fscore_qp <- matrix(0, N_realizations, length(ratios))

print("Connecting to MATLAB...")
matlab <- Matlab()
is_matlab_open <- open(matlab)
print(is_matlab_open)

N <- 64
K <- 4
P <- diag(.49, K) + .01
# K-component graph
A_mask <- matrix(1, N, N) - diag(N)
setVariable(matlab, A_mask = A_mask)

for (j in n_ratios) {
  T <- as.integer(ratios[j] * N)
  cat("\nRunning simulation for", T, "samples per node, p/n = ", ratios[j], "\n")
  for (n in 1:N_realizations) {
    mgraph <- sample_sbm(N, pref.matrix = P, block.sizes = c(rep(N / K, K)))
    E(mgraph)$weight <- runif(gsize(mgraph), min = 1e-1, max = 3)
    Ltrue <- as.matrix(laplacian_matrix(mgraph))
    Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    s_max <- max(abs(S - diag(diag(S))))
    alphas <- c(.75 ^ (c(1:14)) * s_max * sqrt(log(N)/ T), 0)
    setVariable(matlab, S = S)
    Lnaive <- MASS::ginv(S)
    w_qp <- spectralGraphTopology:::w_init("qp", Lnaive)
    Lqp <- L(w_qp)
    graph <- learn_laplacian_matrix(S, w0 = w_qp, ub = 32, beta = 1e2,
                                    maxiter = 1e6, abstol = 0)
    print(graph$convergence)
    rel_cgl <- Inf
    for (alpha in alphas) {
      setVariable(matlab, alpha = alpha)
      evaluate(matlab, "[Lcgl,~,~] = estimate_cgl(S, A_mask, alpha, 1e-4, 1e-4, 100, 1)")
      Lcgl <- getVariable(matlab, "Lcgl")
      if (anyNA(Lcgl$Lcgl))
        next
      tmp_rel_cgl <- relativeError(Ltrue, Lcgl$Lcgl)
      if (tmp_rel_cgl < rel_cgl) {
        rel_cgl <- tmp_rel_cgl
        metrics_cgl <- metrics(Ltrue, Lcgl$Lcgl, eps)
      }
    }
    metrics_sgl <- metrics(Ltrue, graph$Laplacian, eps)
    metrics_qp <- metrics(Ltrue, Lqp, eps)
    metrics_naive <- metrics(Ltrue, Lnaive, eps)
    rel_err_sgl[n, j] <- relativeError(Ltrue, graph$Laplacian)
    rel_err_cgl[n, j] <- rel_cgl
    rel_err_qp[n, j] <- relativeError(Ltrue, Lqp)
    rel_err_naive[n, j] <- relativeError(Ltrue, Lnaive)
    fscore_sgl[n, j] <- metrics_sgl[1]
    fscore_cgl[n, j] <- metrics_cgl[1]
    fscore_qp[n, j] <- metrics_qp[1]
    fscore_naive[n, j] <- metrics_naive[1]
    recall_sgl[n, j] <- metrics_sgl[2]
    recall_cgl[n, j] <- metrics_cgl[2]
    recall_qp[n, j] <- metrics_qp[2]
    recall_naive[n, j] <- metrics_naive[2]
    specificity_sgl[n, j] <- metrics_sgl[3]
    specificity_cgl[n, j] <- metrics_cgl[3]
    specificity_qp[n, j] <- metrics_qp[3]
    specificity_naive[n, j] <- metrics_naive[3]
    accuracy_sgl[n, j] <- metrics_sgl[4]
    accuracy_cgl[n, j] <- metrics_cgl[4]
    accuracy_qp[n, j] <- metrics_qp[4]
    accuracy_naive[n, j] <- metrics_naive[4]
    print(rel_err_sgl)
    print(rel_err_cgl)
    print(fscore_sgl)
    print(fscore_cgl)
  }
}

saveRDS(rel_err_sgl, file = "rel-err-SGL.rds")
saveRDS(rel_err_cgl, file = "rel-err-CGL.rds")
saveRDS(rel_err_naive, file = "rel-err-naive.rds")
saveRDS(rel_err_qp, file = "rel-err-QP.rds")
saveRDS(fscore_sgl, file = "fscore-SGL.rds")
saveRDS(fscore_cgl, file = "fscore-CGL.rds")
saveRDS(fscore_naive, file = "fscore-naive.rds")
saveRDS(fscore_qp, file = "fscore-QP.rds")
saveRDS(recall_sgl, file = "recall-SGL.rds")
saveRDS(recall_cgl, file = "recall-CGL.rds")
saveRDS(recall_naive, file = "recall-naive.rds")
saveRDS(recall_qp, file = "recall-QP.rds")
saveRDS(specificity_sgl, file = "specificity-SGL.rds")
saveRDS(specificity_cgl, file = "specificity-CGL.rds")
saveRDS(specificity_naive, file = "specificity-naive.rds")
saveRDS(specificity_qp, file = "specificity-QP.rds")
saveRDS(accuracy_sgl, file = "accuracy-SGL.rds")
saveRDS(accuracy_cgl, file = "accuracy-CGL.rds")
saveRDS(accuracy_naive, file = "accuracy-naive.rds")
saveRDS(accuracy_qp, file = "accuracy-QP.rds")
