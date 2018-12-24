library(igraph)
library(R.matlab)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)

set.seed(0)
N_realizations <- 1
ratios <- c(2, 10, 100, 1000)
n_ratios <- c(1:length(ratios))
alpha_set <- c(0, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2)
alpha_size <- c(1:length(alpha_set))
alpha_rel_err <- matrix(0, length(alpha_set), length(ratios))
alpha_fscore <- matrix(0, length(alpha_set), length(ratios))

N <- 16
K <- 4
P <- diag(1, K)

for (j in n_ratios) {
  T <- as.integer(ratios[j] * N)
  cat("\nRunning simulation for", T, "samples per node, T/N = ", ratios[j], "\n")
  for (n in 1:N_realizations) {
    mgraph <- sample_sbm(N, pref.matrix = P, block.sizes = c(rep(N / K, K)))
    E(mgraph)$weight <- runif(gsize(mgraph), min = 1e-1, max = 3)
    Ltrue <- as.matrix(laplacian_matrix(mgraph))
    Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    for (i in alpha_size) {
      cat("alpha = ", alpha_set[i], "\n")
      graph <- learnGraphTopology(S, w0 = "qp", K = 4, ub = 2*N,
                                  alpha = alpha_set[i], beta = 10*N, maxiter = 100000)
      print(graph$convergence)
      print(graph$lambda)
      alpha_rel_err[i, j] <- alpha_rel_err[i, j] + relativeError(Ltrue, graph$Lw)
      alpha_fscore[i, j] <- alpha_fscore[i, j] + Fscore(Ltrue, graph$Lw, 1e-1)
    }
  }
  alpha_rel_err[,j] <- alpha_rel_err[,j] / N_realizations
  alpha_fscore[,j] <- alpha_fscore[,j] / N_realizations
}

colors <- c("#0B032D", "#843B62", "#F67E7D", "#FFB997")
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("relative_error_alpha.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(alpha_size, alpha_rel_err[,1], type = "b", pch=15, cex=.75, ylim=c(0, .5),
     xlab = TeX("$\\alpha$"), ylab = "Average Relative Error", col = colors[1], xaxt = "n")
grid()
lines(alpha_size, alpha_rel_err[,2], type = "b", pch=16, cex=.75, col = colors[2], xaxt = "n")
lines(alpha_size, alpha_rel_err[,3], type = "b", pch=17, cex=.75, col = colors[3], xaxt = "n")
lines(alpha_size, alpha_rel_err[,4], type = "b", pch=18, cex=.85, col = colors[4], xaxt = "n")
axis(side = 1, at = alpha_size, labels = alpha_set)
legend("topleft", legend = c(TeX("$T / N = 2$"), TeX("$T / N = 10$"), TeX("$T / N = 100$"),
                              TeX("$T / N = 1000$")),
       col=colors, pch=c(15, 16, 17, 18), lty=c(1, 1, 1, 1), bty="n")
dev.off()
embed_fonts("relative_error_alpha.ps", outfile="relative_error_alpha.ps")
setEPS()
postscript("fscore_alpha.ps", family = "ComputerModern", height = 5, width = gr * 3.5)
plot(alpha_size, alpha_fscore[,1], ylim=c(.5, 1.), xlab = TeX("$\\alpha$"),
     ylab = "Average F-score", type = "b", pch=15, cex=.75, col = colors[1], xaxt = "n")
grid()
lines(alpha_size, alpha_fscore[,2], type = "b", pch=16, cex=.75, col = colors[2], xaxt = "n")
lines(alpha_size, alpha_fscore[,3], type = "b", pch=17, cex=.75, col = colors[3], xaxt = "n")
lines(alpha_size, alpha_fscore[,4], type = "b", pch=18, cex=.85, col = colors[4], xaxt = "n")
axis(side = 1, at = alpha_size, labels = alpha_set)
legend("bottomleft", legend = c(TeX("$T / N = 2$"), TeX("$T / N = 10$"), TeX("$T / N = 100$"),
                                 TeX("$T / N = 1000$")),
       col=colors, pch=c(15, 16, 17, 18), lty=c(1, 1, 1, 1), bty="n")
dev.off()
embed_fonts("fscore_alpha.ps", outfile="fscore_alpha.ps")
