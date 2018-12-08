library(igraph)
library(R.matlab)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)

set.seed(0)
N_realizations <- 5
ratios <- c(1, 10, 100, 1000)
n_ratios <- c(1:length(ratios))
beta_set <- c(.25, 1, 2, 4, 10, 25, 50)
beta_size <- c(1:length(beta_set))
beta_rel_err <- matrix(0, length(beta_set), length(ratios))
beta_fscore <- matrix(0, length(beta_set), length(ratios))

N <- 64
K <- 4
P <- diag(.49, K) + 0.01

for (j in n_ratios) {
  T <- as.integer(ratios[j] * N)
  cat("\nRunning simulation for", T, "samples per node, T/N = ", ratios[j], "\n")
  for (n in 1:N_realizations) {
    mgraph <- sample_sbm(N, pref.matrix = P, block.sizes = c(rep(N / K, K)))
    A_modular_mask <- as.matrix(laplacian_matrix(mgraph)) != 0
    E(mgraph)$weight <- runif(gsize(mgraph), min = 1e-1, max = 3)
    Ltrue <- as.matrix(laplacian_matrix(mgraph))
    Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    for (i in beta_size) {
      print(beta_set[i])
      graph <- learnGraphTopology(S, w0 = "naive", K = 1, ub = 32, beta = beta_set[i], maxiter = 100000)
      print(graph$convergence)
      print(graph$lambda)
      beta_rel_err[i, j] <- beta_rel_err[i, j] + relativeError(Ltrue, graph$Lw)
      beta_fscore[i, j] <- beta_fscore[i, j] + Fscore(Ltrue, graph$Lw, 1e-1)
    }
  }
  beta_rel_err[,j] <- beta_rel_err[,j] / N_realizations
  beta_fscore[,j] <- beta_fscore[,j] / N_realizations
}

colors <- c("#0B032D", "#843B62", "#F67E7D", "#FFB997")
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("relative_error_beta.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(beta_size, beta_rel_err[,1], type = "b", pch=15, cex=.75, ylim=c(0, .5),
     xlab = TeX("$\\beta$"), ylab = "Average Relative Error", col = colors[1], xaxt = "n")
grid()
lines(beta_size, beta_rel_err[,2], type = "b", pch=16, cex=.75, col = colors[2], xaxt = "n")
lines(beta_size, beta_rel_err[,3], type = "b", pch=17, cex=.75, col = colors[3], xaxt = "n")
lines(beta_size, beta_rel_err[,4], type = "b", pch=18, cex=.85, col = colors[4], xaxt = "n")
axis(side = 1, at = beta_size, labels = beta_set)
legend("topleft", legend = c(TeX("$T / N = 1$"), TeX("$T / N = 10$"), TeX("$T / N = 100$"),
                              TeX("$T / N = 1000$")),
       col=colors, pch=c(15, 16, 17, 18), lty=c(1, 1, 1, 1), bty="n")
dev.off()
embed_fonts("relative_error_beta.ps", outfile="relative_error_beta.ps")
setEPS()
postscript("fscore_beta.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(beta_size, beta_fscore[,1], ylim=c(.5, 1.), xlab = TeX("$\\beta$"),
     ylab = "Average F-score", type = "b", pch=15, cex=.75, col = colors[1], xaxt = "n")
grid()
lines(beta_size, beta_fscore[,2], type = "b", pch=16, cex=.75, col = colors[2], xaxt = "n")
lines(beta_size, beta_fscore[,3], type = "b", pch=17, cex=.75, col = colors[3], xaxt = "n")
lines(beta_size, beta_fscore[,4], type = "b", pch=18, cex=.85, col = colors[4], xaxt = "n")
axis(side = 1, at = beta_size, labels = beta_set)
legend("bottomleft", legend = c(TeX("$T / N = 1$"), TeX("$T / N = 10$"), TeX("$T / N = 100$"),
                                 TeX("$T / N = 1000$")),
       col=colors, pch=c(15, 16, 17, 18), lty=c(1, 1, 1, 1), bty="n")
dev.off()
embed_fonts("fscore_beta.ps", outfile="fscore_beta.ps")
