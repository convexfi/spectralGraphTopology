library(igraph)
library(R.matlab)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)
set.seed(23)

N_realizations <- 20
ratios <- c(1., 5., 10, 30, 100, 250, 500, 1000)
n_ratios <- c(1:length(ratios))
# design synthetic Laplacian of a grid graph
N <- 64
grid <- make_lattice(length = sqrt(N), dim = 2)
rel_err_spec <- array(0, length(ratios))
rel_err_cgl <- array(0, length(ratios))
rel_err_naive <- array(0, length(ratios))
rel_err_qp <- array(0, length(ratios))
fscore_spec <- array(0, length(ratios))
fscore_cgl <- array(0, length(ratios))
fscore_naive <- array(0, length(ratios))
fscore_qp <- array(0, length(ratios))

print("Connecting to MATLAB...")
matlab <- Matlab()
open(matlab)
print("success!")
A_mask <- matrix(1, 64, 64) - diag(64)
setVariable(matlab, A_mask = A_mask)

for (j in n_ratios) {
  T <- as.integer(ratios[j] * N)
  cat("\nRunning simulation for", T, "samples per node, T/N = ", ratios[j], "\n")
  for (n in 1:N_realizations) {
    E(grid)$weight <- runif(gsize(grid), min = 1e-1, max = 3)
    Ltrue <- as.matrix(laplacian_matrix(grid))
    # sample data from GP with covariance matrix set as
    # the pseudo inverse of the true Laplacian
    Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    Lnaive <- MASS::ginv(S)
    w_qp <- spectralGraphTopology:::w_init("qp", Lnaive)
    Lqp <- L(w_qp)
    s_max <- max(abs(S - diag(diag(S))))
    alphas <- c(.75 ^ (c(1:14)) * s_max * sqrt(log(N)/ T), 0)
    # run spectralGraphTopology
    if (ratios[j] <= 1)
      graph <- learn_laplacian_matrix(S, w0 = "naive", beta = 1e-2, beta_max = 4, nbeta = 20)#, alpha = 5e-3)
    else
      graph <- learn_laplacian_matrix(S, w0 = w_qp, beta = 10,
                                      alpha = 5e-3, ftol = 1e-6)
    print(graph$convergence)
    print(graph$beta)
    # compute naive
    # set data variable to MATLAB
    setVariable(matlab, S = S)
    rel_cgl <- 9999999999
    for (alpha in alphas) {
      setVariable(matlab, alpha = alpha)
      evaluate(matlab, "[Lcgl,~,~] = estimate_cgl(S, A_mask, alpha, 1e-6, 1e-6, 40, 1)")
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

    rel_spec <- relativeError(Ltrue, graph$Lw)
    fs_spec <- Fscore(Ltrue, graph$Lw, 5e-2)
    print(rel_spec)
    rel_naive <- relativeError(Ltrue, Lnaive)
    fs_naive <- Fscore(Ltrue, Lnaive, 5e-2)
    rel_qp <- relativeError(Ltrue, Lqp)
    fs_qp <- Fscore(Ltrue, Lqp, 5e-2)
    rel_err_spec[j] <- rel_err_spec[j] + rel_spec
    fscore_spec[j] <- fscore_spec[j] + fs_spec
    rel_err_cgl[j] <- rel_err_cgl[j] + rel_cgl
    fscore_cgl[j] <- fscore_cgl[j] + fs_cgl
    rel_err_naive[j] <- rel_err_naive[j] + rel_naive
    fscore_naive[j] <- fscore_naive[j] + fs_naive
    rel_err_qp[j] <- rel_err_qp[j] + rel_qp
    fscore_qp[j] <- fscore_qp[j] + fs_qp

  }
  rel_err_spec[j] <- rel_err_spec[j] / N_realizations
  fscore_spec[j] <- fscore_spec[j] / N_realizations
  rel_err_cgl[j] <- rel_err_cgl[j] / N_realizations
  fscore_cgl[j] <- fscore_cgl[j] / N_realizations
  rel_err_naive[j] <- rel_err_naive[j] / N_realizations
  fscore_naive[j] <- fscore_naive[j] / N_realizations
  rel_err_qp[j] <- rel_err_qp[j] / N_realizations
  fscore_qp[j] <- fscore_qp[j] / N_realizations
  cat("\n** spectralGraphTopology results **\n")
  cat("Avg Relative error: ")
  cat(rel_err_spec[j], "\n")
  cat("Avg Fscore: ")
  cat(fscore_spec[j], "\n")
  cat("\n** CGL results **\n")
  cat("Avg Relative error: ")
  cat(rel_err_cgl[j], "\n")
  cat("Avg Fscore: ")
  cat(fscore_cgl[j], "\n")
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

saveRDS(rel_err_spec, file = "rel-err-SGL.rds")
saveRDS(fscore_spec, file = "fscore-SGL.rds")
saveRDS(rel_err_cgl, file = "rel-err-CGL.rds")
saveRDS(fscore_cgl, file = "fscore-CGL.rds")
saveRDS(rel_err_naive, file = "rel-err-naive.rds")
saveRDS(fscore_naive, file = "fscore-naive.rds")
saveRDS(rel_err_qp, file = "rel-err-QP.rds")
saveRDS(fscore_qp, file = "fscore-QP.rds")

#colors <- c("#0B032D", "#843B62", "#F67E7D", "#FFB997")
#gr = .5 * (1 + sqrt(5))
#setEPS()
#postscript("relative_error_grid.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
#plot(n_ratios, rel_err_naive, type = "b", pch=15, cex=.75, ylim=c(0, .5),
#     xlab = TeX("$\\mathit{T} / \\mathit{N}$"), ylab = "Average Relative Error", col = colors[1], xaxt = "n")
#grid()
#lines(n_ratios, rel_err_cgl, type = "b", pch=16, cex=.75, col = colors[2], xaxt = "n")
#lines(n_ratios, rel_err_spec, type = "b", pch=18, cex=.85, col = colors[4], xaxt = "n")
#axis(side = 1, at = n_ratios, labels = ratios)
#legend("topright", legend = c("ISCM", "CGL", TeX("CGL($\\mathbf{A}$)"), "SGL"),
#       col=colors, pch=c(15, 16, 17, 18), lty=c(1, 1, 1, 1), bty="n")
#dev.off()
#embed_fonts("relative_error_grid.ps", outfile="relative_error_grid.ps")
#setEPS()
#postscript("fscore_grid.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
#plot(n_ratios, fscore_naive, ylim=c(.4, 1.), xlab = TeX("$\\mathit{T} / \\mathit{N}$"),
#     ylab = "Average F-score", type = "b", pch=15, cex=.75, col = colors[1], xaxt = "n")
#grid()
#lines(n_ratios, fscore_cgl, type = "b", pch=16, cex=.75, col = colors[2], xaxt = "n")
#lines(n_ratios, fscore_spec, type = "b", pch=18, cex=.85, col = colors[4], xaxt = "n")
#axis(side = 1, at = n_ratios, labels = ratios)
#legend("bottomright", legend = c("ISCM", "CGL", TeX("CGL($\\mathbf{A}$)"), "SGL"),
#       col=colors, pch=c(15, 16, 17, 18), lty=c(1, 1, 1, 1), bty="n")
#dev.off()
#embed_fonts("fscore_grid.ps", outfile="fscore_grid.ps")
