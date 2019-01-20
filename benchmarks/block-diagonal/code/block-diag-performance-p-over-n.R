library(igraph)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)
library(R.matlab)

set.seed(123)
N_realizations <- 10
ratios <- c(2, 5, 10, 30, 100, 250, 500, 1000)
n_ratios <- c(1:length(ratios))
rel_err_spec <- array(0, length(ratios))
rel_err_naive <- array(0, length(ratios))
rel_err_llqp <- array(0, length(ratios))
fscore_spec <- array(0, length(ratios))
fscore_naive <- array(0, length(ratios))
fscore_llqp <- array(0, length(ratios))

N <- 16
K <- 4
P <- diag(1, K)
# K-component graph
for (j in n_ratios) {
  T <- as.integer(ratios[j] * N)
  cat("\nRunning simulation for", T, "samples per node, p/n = ", ratios[j], "\n")
  for (n in 1:N_realizations) {
    mgraph <- sample_sbm(N, pref.matrix = P, block.sizes = c(rep(N / K, K)))
    E(mgraph)$weight <- runif(gsize(mgraph), min = .1, max = 3)
    A_grid_mask <- as.matrix(laplacian_matrix(mgraph)) != 0
    Ltrue <- as.matrix(laplacian_matrix(mgraph))
    Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    Lnaive <- MASS::ginv(S)
    R <- vecLmat(ncol(Lnaive))
    llqp <- quadprog::solve.QP(crossprod(R), t(R) %*% vec(Lnaive), diag(ncol(R)))
    w0_qp <- llqp$solution
    Lqp <- L(w0_qp)
    s_max <- max(abs(S - diag(diag(S))))
    alphas <- c(0, .75 ^ (c(1:14)) * s_max * sqrt(log(N)/ T))
    rel_spec <- 9999999999
    for (alpha in alphas) {
      graph <- learnLaplacianGraphTopology(S, w0 = w0_qp, K = K, beta = 10*N,
                                           alpha = alpha, ub = 2*N, maxiter = 1e5)
      print(graph$convergence)
      tmp_rel_spec <- relativeError(Ltrue, graph$Lw)
      if (tmp_rel_spec < rel_spec) {
        fs_spec <- Fscore(Ltrue, graph$Lw, 1e-1)
        rel_spec = tmp_rel_spec
        cat("\nalpha = ", alpha)
      }
    }
    rel_naive = relativeError(Ltrue, Lnaive)
    fs_naive = Fscore(Ltrue, Lnaive, 1e-1)
    rel_llqp <- relativeError(Ltrue, Lqp)
    fs_llqp <- Fscore(Ltrue, Lqp, 1e-1)
    rel_err_spec[j] <- rel_err_spec[j] + rel_spec
    fscore_spec[j] <- fscore_spec[j] + fs_spec
    rel_err_naive[j] <- rel_err_naive[j] + rel_naive
    fscore_naive[j] <- fscore_naive[j] + fs_naive
    rel_err_llqp[j] <- rel_err_llqp[j] + rel_llqp
    fscore_llqp[j] <- fscore_llqp[j] + fs_llqp
  }
  rel_err_spec[j] <- rel_err_spec[j] / N_realizations
  fscore_spec[j] <- fscore_spec[j] / N_realizations
  rel_err_naive[j] <- rel_err_naive[j] / N_realizations
  fscore_naive[j] <- fscore_naive[j] / N_realizations
  rel_err_llqp[j] <- rel_err_llqp[j] / N_realizations
  fscore_llqp[j] <- fscore_llqp[j] / N_realizations
  cat("\n** spectralGraphTopology results **\n")
  cat("Avg Relative error: ")
  cat(rel_err_spec[j], "\n")
  cat("Avg Fscore: ")
  cat(fscore_spec[j], "\n")
  cat("\n** Naive results **\n")
  cat("Avg Relative error: ")
  cat(rel_err_naive[j], "\n")
  cat("Avg Fscore: ")
  cat(fscore_naive[j], "\n")
  cat("\n** LLQP results **\n")
  cat("Avg Relative error: ")
  cat(rel_err_llqp[j], "\n")
  cat("Avg Fscore: ")
  cat(fscore_llqp[j], "\n")
}

saveRDS(rel_err_spec, file="rel-err-SGL.rds")
saveRDS(rel_err_naive, file="rel-err-Naive.rds")
saveRDS(rel_err_llqp, file="rel-err-QP.rds")
saveRDS(fscore_spec, file="fscore-SGL.rds")
saveRDS(fscore_naive, file="fscore-Naive.rds")
saveRDS(fscore_llqp, file="fscore-QP.rds")
