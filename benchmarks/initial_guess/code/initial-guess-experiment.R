library(igraph)
library(R.matlab)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)

N_realizations <- 10
ratios <- c(1, 3, 5, 10, 20, 50, 100)
n_ratios <- c(1:length(ratios))
# design synthetic Laplacian of a grid graph
N <- 64
grid <- make_lattice(length = sqrt(N), dim = 2)
rel_err_spec_naive <- array(0, length(ratios))
rel_err_naive <- array(0, length(ratios))
rel_err_spec_qp <- array(0, length(ratios))
rel_err_qp <- array(0, length(ratios))
fscore_spec_naive <- array(0, length(ratios))
fscore_naive <- array(0, length(ratios))
fscore_spec_qp <- array(0, length(ratios))
fscore_qp <- array(0, length(ratios))

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
    R <- vecLmat(ncol(Lnaive))
    qp <- quadprog::solve.QP(crossprod(R), t(R) %*% vec(Lnaive), diag(ncol(R)))
    w0_qp <- qp$solution
    Lqp <- L(w0_qp)
    # run spectralGraphTopology
    if (ratios[j] <= 5) {
      graph <- learnLaplacianGraphTopology(S, w0 = "naive", beta = 1e-2, beta_max = 4, nbeta = 20)
      graph_qp <- learnLaplacianGraphTopology(S, w0 = w0_qp, beta = 1e-2, beta_max = 4, nbeta = 20)
    } else {
      graph <- learnLaplacianGraphTopology(S, w0 = "naive", beta = 10)
      graph_qp <- learnLaplacianGraphTopology(S, w0 = w0_qp, beta = 10)
    }
    rel_spec_naive <- relativeError(Ltrue, graph$Lw)
    fs_spec_naive <- Fscore(Ltrue, graph$Lw, 1e-1)
    rel_spec_qp <- relativeError(Ltrue, graph_qp$Lw)
    fs_spec_qp <- Fscore(Ltrue, graph_qp$Lw, 1e-1)
    rel_qp <- relativeError(Ltrue, Lqp)
    fs_qp <- Fscore(Ltrue, Lqp, 1e-1)
    rel_naive <- relativeError(Ltrue, Lnaive)
    fs_naive <- Fscore(Ltrue, Lnaive, 1e-1)

    rel_err_spec_naive[j] <- rel_err_spec_naive[j] + rel_spec_naive
    fscore_spec_naive[j] <- fscore_spec_naive[j] + fs_spec_naive
    rel_err_spec_qp[j] <- rel_err_spec_qp[j] + rel_spec_qp
    fscore_spec_qp[j] <- fscore_spec_qp[j] + fs_spec_qp

    rel_err_naive[j] <- rel_err_naive[j] + rel_naive
    fscore_naive[j] <- fscore_naive[j] + fs_naive
    rel_err_qp[j] <- rel_err_qp[j] + rel_qp
    fscore_qp[j] <- fscore_qp[j] + fs_qp
  }
  rel_err_spec_naive[j] <- rel_err_spec_naive[j] / N_realizations
  fscore_spec_naive[j] <- fscore_spec_naive[j] / N_realizations
  rel_err_spec_qp[j] <- rel_err_spec_qp[j] / N_realizations
  fscore_spec_qp[j] <- fscore_spec_qp[j] / N_realizations
  rel_err_naive[j] <- rel_err_naive[j] / N_realizations
  fscore_naive[j] <- fscore_naive[j] / N_realizations
  rel_err_qp[j] <- rel_err_qp[j] / N_realizations
  fscore_qp[j] <- fscore_qp[j] / N_realizations
  cat("\n** Spec with Naive results **\n")
  cat("Avg Relative error: ")
  cat(rel_err_spec_naive[j], "\n")
  cat("Avg Fscore: ")
  cat(fscore_spec_naive[j], "\n")
  cat("\n** Spec with QP results **\n")
  cat("Avg Relative error: ")
  cat(rel_err_spec_qp[j], "\n")
  cat("Avg Fscore: ")
  cat(fscore_spec_qp[j], "\n")
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

colors <- c("#706FD3", "#FF5252", "#33D9B2", "#34ACE0")
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("../latex/figures/relative_error.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(n_ratios, rel_err_spec_naive, type = "b", pch=15, cex=.75, ylim=c(0, .7),
     xlab = TeX("$\\mathit{T} / \\mathit{N}$"), ylab = "Average Relative Error", col = colors[1], xaxt = "n")
grid()
lines(n_ratios, rel_err_spec_qp, type = "b", pch=16, cex=.75, col = colors[2], xaxt = "n")
lines(n_ratios, rel_err_naive, type = "b", pch=17, cex=.75, col = colors[3], xaxt = "n")
lines(n_ratios, rel_err_qp, type = "b", pch=18, cex=.85, col = colors[4], xaxt = "n")
axis(side = 1, at = n_ratios, labels = ratios)
legend("topright", legend = c(TeX("SGL($w_0 = 'naive'$)"), TeX("SGL($w_0 = 'qp'$)"), "ISCM", "LLQP"),
       col=colors, pch=c(15, 16, 17, 18), lty=c(1, 1, 1, 1), bty="n")
dev.off()
embed_fonts("../latex/figures/relative_error.ps", outfile="../latex/figures/relative_error.ps")
setEPS()
postscript("../latex/figures/fscore.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(n_ratios, fscore_spec_naive, ylim=c(.1, 1.), xlab = TeX("$\\mathit{T} / \\mathit{N}$"),
     ylab = "Average F-score", type = "b", pch=15, cex=.75, col = colors[1], xaxt = "n")
grid()
lines(n_ratios, fscore_spec_qp, type = "b", pch=16, cex=.75, col = colors[2], xaxt = "n")
lines(n_ratios, fscore_naive, type = "b", pch=17, cex=.75, col = colors[3], xaxt = "n")
lines(n_ratios, fscore_qp, type = "b", pch=18, cex=.85, col = colors[4], xaxt = "n")
axis(side = 1, at = n_ratios, labels = ratios)
legend("bottomright", legend = c(TeX("SGL($w_0 = 'naive'$)"), TeX("SGL($w_0 = 'qp'$)"), "ISCM", "LLQP"),
       col=colors, pch=c(15, 16, 17, 18), lty=c(1, 1, 1, 1), bty="o", bg = "white", box.col = "white")
dev.off()
embed_fonts("../latex/figures/fscore.ps", outfile="../latex/figures/fscore.ps")

