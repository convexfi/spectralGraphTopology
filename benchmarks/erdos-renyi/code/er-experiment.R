library(igraph)
library(R.matlab)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)

set.seed(42)
eps <- 5e-2
N_realizations <- 2
ratios <- rev(c(.5, 5, 10, 50, 100, 500, 1e3, 5e3))
n_ratios <- c(1:length(ratios))
# design synthetic Laplacian of a erdos_renyi graph
N <- 64
p <- .25
rel_err_sgl <- matrix(0, N_realizations, length(ratios))
rel_err_cgl <- matrix(0, N_realizations, length(ratios))
rel_err_qp <- matrix(0, N_realizations, length(ratios))
rel_err_naive <- matrix(0, N_realizations, length(ratios))
fscore_sgl <- matrix(0, N_realizations, length(ratios))
fscore_cgl <- matrix(0, N_realizations, length(ratios))
fscore_qp <- matrix(0, N_realizations, length(ratios))
fscore_naive <- matrix(0, N_realizations, length(ratios))
accuracy_sgl <- matrix(0, N_realizations, length(ratios))
accuracy_cgl <- matrix(0, N_realizations, length(ratios))
accuracy_qp <- matrix(0, N_realizations, length(ratios))
accuracy_naive <- matrix(0, N_realizations, length(ratios))
specificity_sgl <- matrix(0, N_realizations, length(ratios))
specificity_cgl <- matrix(0, N_realizations, length(ratios))
specificity_qp <- matrix(0, N_realizations, length(ratios))
specificity_naive <- matrix(0, N_realizations, length(ratios))

print("Connecting to MATLAB...")
matlab <- Matlab(port=9998)
open(matlab)
print("success!")
A_mask <- matrix(1, N, N) - diag(N)
setVariable(matlab, A_mask = A_mask)

for (j in n_ratios) {
  T <- as.integer(ratios[j] * N)
  cat("\nRunning simulation for", T, "samples per node, T/N = ", ratios[j], "\n")
  for (n in 1:N_realizations) {
    erdos_renyi <- erdos.renyi.game(N, p)
    A_erdos_renyi_mask <- as.matrix(laplacian_matrix(erdos_renyi)) != 0
    setVariable(matlab, A_erdos_renyi_mask = A_erdos_renyi_mask)
    E(erdos_renyi)$weight <- runif(gsize(erdos_renyi), min = 1e-1, max = 3)
    Ltrue <- as.matrix(laplacian_matrix(erdos_renyi))
    # sample data from GP with covariance matrix set as
    # the pseudo inverse of the true Laplacian
    Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    Sinv <- MASS::ginv(S)
    w_qp <- spectralGraphTopology:::w_init("qp", Sinv)
    w_naive <- spectralGraphTopology:::w_init("naive", Sinv)
    Lnaive <- L(w_naive)
    Lqp <- L(w_qp)
    s_max <- max(abs(S - diag(diag(S))))
    alphas <- c(.75 ^ (c(1:14)) * s_max * sqrt(log(N)/T), 0)
    setVariable(matlab, S = S)
    rel_cgl <- Inf
    for (alpha in alphas) {
      setVariable(matlab, alpha = alpha)
      evaluate(matlab, "[Lcgl,~,~] = estimate_cgl(S, A_mask, alpha, 1e-4, 1e-4, 100, 1)")
      Lcgl <- getVariable(matlab, "Lcgl")
      if (anyNA(Lcgl$Lcgl)) {
        next
      }
      tmp_rel_cgl <- relativeError(Ltrue, Lcgl$Lcgl)
      if (tmp_rel_cgl < rel_cgl) {
        alpha_tmp <- alpha
        rel_cgl <- tmp_rel_cgl
        metrics_cgl <- metrics(Ltrue, Lcgl$Lcgl, eps)
      }
    }
    print(alpha_tmp)
    graph <- learn_k_component_graph(S, w0 = w_qp, beta = 1, ub = 32, fix_beta = TRUE,
                                     alpha = 0, maxiter = 5e5, abstol = 0)
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
saveRDS(specificity_sgl, file = "specificity-SGL.rds")
saveRDS(specificity_cgl, file = "specificity-CGL.rds")
saveRDS(specificity_naive, file = "specificity-naive.rds")
saveRDS(specificity_qp, file = "specificity-QP.rds")
saveRDS(accuracy_sgl, file = "accuracy-SGL.rds")
saveRDS(accuracy_cgl, file = "accuracy-CGL.rds")
saveRDS(accuracy_naive, file = "accuracy-naive.rds")
saveRDS(accuracy_qp, file = "accuracy-QP.rds")
