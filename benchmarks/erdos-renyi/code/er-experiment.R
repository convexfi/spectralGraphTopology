library(igraph)
library(R.matlab)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)

set.seed(42)

N_realizations <- 2
ratios <- c(.5, .75, 1.5, 5, 10, 30, 50, 100, 250, 500, 1000)
n_ratios <- c(1:length(ratios))
# design synthetic Laplacian of a erdos_renyi graph
N <- 64
p <- .4
rel_err_spec <- matrix(0, N_realizations, length(ratios))
rel_err_cgl <- matrix(0, N_realizations, length(ratios))
rel_err_qp <- matrix(0, N_realizations, length(ratios))
rel_err_naive <- matrix(0, N_realizations, length(ratios))
fscore_spec <- matrix(0, N_realizations, length(ratios))
fscore_cgl <- matrix(0, N_realizations, length(ratios))
fscore_qp <- matrix(0, N_realizations, length(ratios))
fscore_naive <- matrix(0, N_realizations, length(ratios))

print("Connecting to MATLAB...")
matlab <- Matlab(port=9999)
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
    Lnaive <- MASS::ginv(S)
    w_qp <- spectralGraphTopology:::w_init("qp", Lnaive)
    Lqp <- L(w_qp)
    s_max <- max(abs(S - diag(diag(S))))
    alphas <- c(.75 ^ (c(1:14)) * s_max * sqrt(log(N)/T), 0)
    if (ratios[j] <= 10)
      graph <- learn_laplacian_matrix(S, w0 = w_qp, beta = 1, tol = 1e-6, maxiter = 1e6)
    else
      graph <- learn_laplacian_matrix(S, w0 = w_qp, beta = 100, tol = 1e-6, maxiter = 1e6)
    print(graph$beta_seq)
    print(graph$convergence)
    setVariable(matlab, S = S)
    rel_cgl <- Inf
    for (alpha in alphas) {
      setVariable(matlab, alpha = alpha)
      evaluate(matlab, "[Lcgl,~,~] = estimate_cgl(S, A_mask, alpha, 1e-6, 1e-6, 100, 1)")
      Lcgl <- getVariable(matlab, "Lcgl")
      if (anyNA(Lcgl$Lcgl)) {
        next
      }
      tmp_rel_cgl <- relativeError(Ltrue, Lcgl$Lcgl)
      if (tmp_rel_cgl < rel_cgl) {
        rel_cgl <- tmp_rel_cgl
        fs_cgl <- Fscore(Ltrue, Lcgl$Lcgl, 5e-2)
      }
    }
    rel_err_spec[n, j] <- relativeError(Ltrue, graph$Laplacian)
    fscore_spec[n, j] <- Fscore(Ltrue, graph$Laplacian, 5e-2)
    print(rel_err_spec)
    print(fscore_spec)
    rel_err_cgl[n, j] <- rel_cgl
    fscore_cgl[n, j] <- fs_cgl
    print(rel_err_cgl)
    print(fscore_cgl)
    rel_err_qp[n, j] <- relativeError(Ltrue, Lqp)
    fscore_qp[n, j] <- Fscore(Ltrue, Lqp, 5e-2)
    rel_err_naive[n, j] <- relativeError(Ltrue, Lnaive)
    fscore_naive[n, j] <- Fscore(Ltrue, Lnaive, 5e-2)
  }
}

saveRDS(rel_err_spec, file = "rel-err-SGL.rds")
saveRDS(fscore_spec, file = "fscore-SGL.rds")
saveRDS(rel_err_cgl, file = "rel-err-CGL.rds")
saveRDS(fscore_cgl, file = "fscore-CGL.rds")
saveRDS(rel_err_naive, file = "rel-err-naive.rds")
saveRDS(fscore_naive, file = "fscore-naive.rds")
saveRDS(rel_err_qp, file = "rel-err-QP.rds")
saveRDS(fscore_qp, file = "fscore-QP.rds")
