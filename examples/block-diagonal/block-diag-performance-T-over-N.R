library(igraph)
library(R.matlab)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)

N_realizations <- 1
ratios <- c(500, 1000)
n_ratios <- c(1:length(ratios))
rel_err_spec <- array(0, length(ratios))
rel_err_cglA <- array(0, length(ratios))
rel_err_naive <- array(0, length(ratios))
fscore_spec <- array(0, length(ratios))
fscore_cglA <- array(0, length(ratios))
fscore_naive <- array(0, length(ratios))

print("Connecting to MATLAB...")
matlab <- Matlab(port=9997)
is_matlab_open <- open(matlab)
print(is_matlab_open)

N <- 64
K <- 4
P <- diag(1, K)
# K-component graph
mgraph <- sample_sbm(N, pref.matrix = P, block.sizes = c(rep(N / K, K)))
A_mask <- as.matrix(laplacian_matrix(mgraph)) != 0
setVariable(matlab, A_mask = A_mask)

for (j in n_ratios) {
  fails <- 0
  T <- as.integer(ratios[j] * N)
  cat("\nRunning simulation for", T, "samples per node, T/N = ", ratios[j], "\n")
  for (n in 1:N_realizations) {
    E(mgraph)$weight <- runif(gsize(mgraph), min = 1e-1, max = 3)
    Ltrue <- as.matrix(laplacian_matrix(mgraph))
    Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y) + (1/3) * diag(N)
    s_max <- max(abs(S - diag(diag(S))))
    alphas <- c(.75 ^ (c(1:14)) * s_max * sqrt(log(N)/ T), 0)
    setVariable(matlab, S = S)
    if (ratios[j] <= 5) {
      graph <- learnGraphTopology(S, w0 = "naive", K = K, beta = 1e-3, beta_max = 2, nbeta = 20)
    } else {
      graph <- learnGraphTopology(S, w0 = "naive", K = K, beta = 2)
    }
    rel_cgl <- 9999999999
    rel_cglA <- 9999999999
    for (alpha in alphas) {
      setVariable(matlab, alpha = alpha)
      evaluate(matlab, "[LcglA,~,~] = estimate_cgl(S, A_mask, alpha, 1e-4, 1e-6, 40, 1)")
      LcglA <- getVariable(matlab, "LcglA")
      if (anyNA(LcglA$LcglA)) {
        next
      }
      tmp_rel_cglA <- relativeError(Ltrue, LcglA$LcglA)
      if (tmp_rel_cglA < rel_cglA) {
        rel_cglA <- tmp_rel_cglA
        fs_cglA <- Fscore(Ltrue, LcglA$LcglA, 1e-1)
      }
    }
    Lnaive <- MASS::ginv(S)
    # set data variable to MATLAB
    rel_spec = relativeError(Ltrue, graph$Lw)
    fs_spec = Fscore(Ltrue, graph$Lw, 1e-1)
    rel_naive = relativeError(Ltrue, Lnaive)
    fs_naive = Fscore(Ltrue, Lnaive, 1e-1)
    print(rel_spec)
    print(fs_spec)
    print(rel_cglA)
    print(fs_cglA)
    print(rel_naive)
    print(fs_naive)
    rel_err_spec[j] <- rel_err_spec[j] + rel_spec
    fscore_spec[j] <- fscore_spec[j] + fs_spec
    rel_err_cglA[j] <- rel_err_cglA[j] + rel_cglA
    fscore_cglA[j] <- fscore_cglA[j] + fs_cglA
    rel_err_naive[j] <- rel_err_naive[j] + rel_naive
    fscore_naive[j] <- fscore_naive[j] + fs_naive
  }
  rel_err_spec[j] <- rel_err_spec[j] / N_realizations
  fscore_spec[j] <- fscore_spec[j] / N_realizations
  rel_err_cglA[j] <- rel_err_cglA[j] / (N_realizations - fails)
  fscore_cglA[j] <- fscore_cglA[j] / (N_realizations - fails)
  rel_err_naive[j] <- rel_err_naive[j] / N_realizations
  fscore_naive[j] <- fscore_naive[j] / N_realizations
  cat("\n** spectralGraphTopology results **\n")
  cat("Avg Relative error: ")
  cat(rel_err_spec[j], "\n")
  cat("Avg Fscore: ")
  cat(fscore_spec[j], "\n")
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

colors <- c("#0B032D", "#F67E7D", "#FFB997")
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("relative_error_block_diagonal.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(n_ratios, rel_err_naive, type = "b", pch=15, cex=.75, ylim=c(0, .5),
     xlab = TeX("$\\mathit{T} / \\mathit{N}$"), ylab = "Average Relative Error", col = colors[1], xaxt = "n")
grid()
lines(n_ratios, rel_err_cglA, type = "b", pch=17, cex=.75, col = colors[2], xaxt = "n")
lines(n_ratios, rel_err_spec, type = "b", pch=18, cex=.85, col = colors[3], xaxt = "n")
axis(side = 1, at = n_ratios, labels = ratios)
legend("topright", legend = c("ISCM", TeX("CGL($\\mathbf{A}$)"), "SGL"),
       col=colors, pch=c(15, 17, 18), lty=c(1, 1, 1), bty="n")
dev.off()
embed_fonts("relative_error_block_diagonal.ps", outfile="relative_error_block_diagonal.ps")
setEPS()
postscript("fscore_block_diagonal.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(n_ratios, fscore_naive, ylim=c(.3, 1.), xlab = TeX("$\\mathit{T} / \\mathit{N}$"),
     ylab = "Average F-score", type = "b", pch=15, cex=.75, col = colors[1], xaxt = "n")
grid()
lines(n_ratios, fscore_cglA, type = "b", pch=17, cex=.75, col = colors[2], xaxt = "n")
lines(n_ratios, fscore_spec, type = "b", pch=18, cex=.85, col = colors[3], xaxt = "n")
axis(side = 1, at = n_ratios, labels = ratios)
legend("bottomright", legend = c("ISCM", TeX("CGL($\\mathbf{A}$)"), "SGL"),
       col=colors, pch=c(15, 17, 18), lty=c(1, 1, 1), bty="n")
dev.off()
embed_fonts("fscore_block_diagonal.ps", outfile="fscore_block_diagonal.ps")
