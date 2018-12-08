library(igraph)
library(R.matlab)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)

set.seed(12)
N_realizations <- 20
ratios <- c(.5, .75, 1., 5., 10, 30, 100, 250, 500, 1000)
n_ratios <- c(1:length(ratios))
rel_err_spec <- array(0, length(ratios))
rel_err_cgl <- array(0, length(ratios))
rel_err_cglA <- array(0, length(ratios))
rel_err_naive <- array(0, length(ratios))
fscore_spec <- array(0, length(ratios))
fscore_cgl <- array(0, length(ratios))
fscore_cglA <- array(0, length(ratios))
fscore_naive <- array(0, length(ratios))

print("Connecting to MATLAB...")
matlab <- Matlab(port=9999)
is_matlab_open <- open(matlab)
print(is_matlab_open)

N <- 64
K <- 4
P <- diag(.49, K) + 0.01
# K-component graph
A_mask <- matrix(1, N, N) - diag(N)
setVariable(matlab, A_mask = A_mask)

for (j in n_ratios) {
  T <- as.integer(ratios[j] * N)
  cat("\nRunning simulation for", T, "samples per node, T/N = ", ratios[j], "\n")
  for (n in 1:N_realizations) {
    mgraph <- sample_sbm(N, pref.matrix = P, block.sizes = c(rep(N / K, K)))
    A_modular_mask <- as.matrix(laplacian_matrix(mgraph)) != 0
    setVariable(matlab, A_modular_mask = A_modular_mask)
    E(mgraph)$weight <- runif(gsize(mgraph), min = 1e-1, max = 3)
    Ltrue <- as.matrix(laplacian_matrix(mgraph))
    Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    s_max <- max(abs(S - diag(diag(S))))
    alphas <- c(.75 ^ (c(1:14)) * s_max * sqrt(log(N)/ T), 0)
    setVariable(matlab, S = S)
    graph <- learnGraphTopology(S, w0 = "naive", K = 1, ub = 32, beta = 4, maxiter = 100000)
    print(graph$convergence)
    print(graph$lambda)
    rel_cgl <- 9999999999
    rel_cglA <- 9999999999
    for (alpha in alphas) {
      setVariable(matlab, alpha = alpha)
      evaluate(matlab, "[Lcgl,~,~] = estimate_cgl(S, A_mask, alpha, 1e-6, 1e-6, 40, 1)")
      evaluate(matlab, "[LcglA,~,~] = estimate_cgl(S, A_modular_mask, alpha, 1e-6, 1e-6, 40, 1)")
      Lcgl <- getVariable(matlab, "Lcgl")
      LcglA <- getVariable(matlab, "LcglA")
      if (anyNA(Lcgl$Lcgl) || anyNA(LcglA$LcglA)) {
        next
      }
      tmp_rel_cgl <- relativeError(Ltrue, Lcgl$Lcgl)
      if (tmp_rel_cgl < rel_cgl) {
        rel_cgl <- tmp_rel_cgl
        fs_cgl <- Fscore(Ltrue, Lcgl$Lcgl, 1e-1)
      }
      tmp_rel_cglA <- relativeError(Ltrue, LcglA$LcglA)
      if (tmp_rel_cglA < rel_cglA) {
        rel_cglA <- tmp_rel_cglA
        fs_cglA <- Fscore(Ltrue, LcglA$LcglA, 1e-1)
      }
    }
    Lnaive <- MASS::ginv(S)
    rel_spec = relativeError(Ltrue, graph$Lw)
    fs_spec = Fscore(Ltrue, graph$Lw, 1e-1)
    rel_naive = relativeError(Ltrue, Lnaive)
    fs_naive = Fscore(Ltrue, Lnaive, 1e-1)
    rel_err_spec[j] <- rel_err_spec[j] + rel_spec
    fscore_spec[j] <- fscore_spec[j] + fs_spec
    rel_err_cglA[j] <- rel_err_cglA[j] + rel_cglA
    fscore_cglA[j] <- fscore_cglA[j] + fs_cglA
    rel_err_cgl[j] <- rel_err_cgl[j] + rel_cgl
    fscore_cgl[j] <- fscore_cgl[j] + fs_cgl
    rel_err_naive[j] <- rel_err_naive[j] + rel_naive
    fscore_naive[j] <- fscore_naive[j] + fs_naive
  }
  rel_err_spec[j] <- rel_err_spec[j] / N_realizations
  fscore_spec[j] <- fscore_spec[j] / N_realizations
  rel_err_cglA[j] <- rel_err_cglA[j] / N_realizations
  fscore_cglA[j] <- fscore_cglA[j] / N_realizations
  rel_err_cgl[j] <- rel_err_cgl[j] / N_realizations
  fscore_cgl[j] <- fscore_cgl[j] / N_realizations
  rel_err_naive[j] <- rel_err_naive[j] / N_realizations
  fscore_naive[j] <- fscore_naive[j] / N_realizations
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
  cat("\n** CGL with A results **\n")
  cat("Avg Relative error: ")
  cat(rel_err_cglA[j], "\n")
  cat("Avg Fscore: ")
  cat(fscore_cglA[j], "\n")
  cat("\n** Naive results **\n")
  cat("Avg Relative error: ")
  cat(rel_err_naive[j], "\n")
  cat("Avg Fscore: ")
  cat(fscore_naive[j], "\n")
}

colors <- c("#0B032D", "#843B62", "#F67E7D", "#FFB997")
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("relative_error_modular.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(n_ratios, rel_err_naive, type = "b", pch=15, cex=.75, ylim=c(0, .5),
     xlab = TeX("$\\mathit{T} / \\mathit{N}$"), ylab = "Average Relative Error", col = colors[1], xaxt = "n")
grid()
lines(n_ratios, rel_err_cgl, type = "b", pch=16, cex=.75, col = colors[2], xaxt = "n")
lines(n_ratios, rel_err_cglA, type = "b", pch=17, cex=.75, col = colors[3], xaxt = "n")
lines(n_ratios, rel_err_spec, type = "b", pch=18, cex=.85, col = colors[4], xaxt = "n")
axis(side = 1, at = n_ratios, labels = ratios)
legend("topright", legend = c("ISCM", "CGL", TeX("CGL($\\mathbf{A}$)"), "SGL"),
       col=colors, pch=c(15, 16, 17, 18), lty=c(1, 1, 1, 1), bty="n")
dev.off()
embed_fonts("relative_error_modular.ps", outfile="relative_error_modular.ps")
setEPS()
postscript("fscore_modular.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(n_ratios, fscore_naive, ylim=c(.5, 1.), xlab = TeX("$\\mathit{T} / \\mathit{N}$"),
     ylab = "Average F-score", type = "b", pch=15, cex=.75, col = colors[1], xaxt = "n")
grid()
lines(n_ratios, fscore_cgl, type = "b", pch=16, cex=.75, col = colors[2], xaxt = "n")
lines(n_ratios, fscore_cglA, type = "b", pch=17, cex=.75, col = colors[3], xaxt = "n")
lines(n_ratios, fscore_spec, type = "b", pch=18, cex=.85, col = colors[4], xaxt = "n")
axis(side = 1, at = n_ratios, labels = ratios)
legend("bottomright", legend = c("ISCM", "CGL", TeX("CGL($\\mathbf{A}$)"), "SGL"),
       col=colors, pch=c(15, 16, 17, 18), lty=c(1, 1, 1, 1), bty="o", bg = "white", box.col = "white")
dev.off()
embed_fonts("fscore_modular.ps", outfile="fscore_modular.ps")
