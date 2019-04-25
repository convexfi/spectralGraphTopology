library(igraph)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)

set.seed(42)
N_realizations <- 10
ratios <- c(10, 100, 1000)
n_ratios <- c(1:length(ratios))
beta_set <- c(1, 10, 20, 50, 1e2, 1e3, 1e4, 1e5)
beta_size <- c(1:length(beta_set))
beta_rel_err <- matrix(0, length(beta_set), length(ratios))
beta_fscore <- matrix(0, length(beta_set), length(ratios))

N <- 32
K <- 4
P <- diag(.5, K)

mgraph <- sample_sbm(N, pref.matrix = P, block.sizes = c(rep(N / K, K)))
for (j in n_ratios) {
  T <- as.integer(ratios[j] * N)
  cat("\nRunning simulation for", T, "samples per node, n / p = ", ratios[j], "\n")
  for (n in 1:N_realizations) {
    E(mgraph)$weight <- runif(gsize(mgraph), min = 0, max = 1)
    Ltrue <- as.matrix(laplacian_matrix(mgraph))
    Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    for (i in beta_size) {
      print(beta_set[i])
      graph <- learn_k_component_graph(S, w0 = "qp", k = 4, ub = 32, beta = beta_set[i],
                                      fix_beta = TRUE, maxiter = 100000)
      beta_rel_err[i, j] <- beta_rel_err[i, j] + relativeError(Ltrue, graph$Laplacian)
      beta_fscore[i, j] <- beta_fscore[i, j] + metrics(Ltrue, graph$Laplacian, 5e-2)[1]
    }
  }
  beta_rel_err[,j] <- beta_rel_err[,j] / N_realizations
  beta_fscore[,j] <- beta_fscore[,j] / N_realizations
}

saveRDS(beta_rel_err, file = "beta_rel_err.rds")
saveRDS(beta_fscore, file = "beta_fscore.rds")
#colors <- c("#0B032D", "#843B62", "#F67E7D", "#FFB997")
#gr = .5 * (1 + sqrt(5))
#setEPS()
#postscript("relative_error_beta.ps", family = "ComputerModern", height = 5.5, width = gr * 4.5, pointsize = 14)
#plot(beta_size, beta_rel_err[,1], type = "b", pch=15, cex=.75, ylim=c(0, max(beta_rel_err)),
#     xlab = TeX("$\\beta$"), ylab = "Average Relative Error", col = colors[1], xaxt = "n")
#grid()
#lines(beta_size, beta_rel_err[,2], type = "b", pch=16, cex=.75, col = colors[2], xaxt = "n")
#lines(beta_size, beta_rel_err[,3], type = "b", pch=17, cex=.75, col = colors[3], xaxt = "n")
#lines(beta_size, beta_rel_err[,4], type = "b", pch=18, cex=.85, col = colors[4], xaxt = "n")
#axis(side = 1, at = beta_size, labels = beta_set)
#legend("topright", legend = c(TeX("$n / p = .5$"), TeX("$n / p = 10$"), TeX("$n / p = 100$"), TeX("$n / p = 1000$")),
#       col=colors, pch=c(15, 16, 17, 18), lty=c(1, 1, 1, 1), bty="n")
#dev.off()
#embed_fonts("relative_error_beta.ps", outfile="relative_error_beta.ps")
#setEPS()
#postscript("fscore_beta.ps", family = "ComputerModern", height = 5.5, width = gr * 4.5, pointsize = 14)
#plot(beta_size, beta_fscore[,1], ylim=c(min(beta_fscore), 1.), xlab = TeX("$\\beta$"),
#     ylab = "Average F-score", type = "b", pch=15, cex=.75, col = colors[1], xaxt = "n")
#grid()
#lines(beta_size, beta_fscore[,2], type = "b", pch=16, cex=.75, col = colors[2], xaxt = "n")
#lines(beta_size, beta_fscore[,3], type = "b", pch=17, cex=.75, col = colors[3], xaxt = "n")
#lines(beta_size, beta_fscore[,4], type = "b", pch=18, cex=.85, col = colors[4], xaxt = "n")
#axis(side = 1, at = beta_size, labels = beta_set)
#legend("topleft", legend = c(TeX("$n / p = .5$"), TeX("$n / p = 10$"), TeX("$n / p = 100$"), TeX("$n / p = 1000$")),
#       col=colors, pch=c(15, 16, 17, 18), lty=c(1, 1, 1, 1), bty="n")
#dev.off()
#embed_fonts("fscore_beta.ps", outfile="fscore_beta.ps")
