library(igraph)
library(spectralGraphTopology)
library(extrafont)
library(latex2exp)

set.seed(42)
N_realizations <- 10
ratios <- c(10, 100, 1000)
n_ratios <- c(1:length(ratios))
gamma_set <- c(1, 10, 20, 50, 1e2, 1e3, 1e4, 1e5)
gamma_size <- c(1:length(gamma_set))
gamma_rel_err <- matrix(0, length(gamma_set), length(ratios))
gamma_fscore <- matrix(0, length(gamma_set), length(ratios))

n1 <- 20
n2 <- 12
b1 <- sample_bipartite(n1, n2, type="Gnp", p = .7, directed=FALSE)
for (j in n_ratios) {
  T <- as.integer(ratios[j] * (n1 + n2))
  cat("\nRunning simulation for", T, "samples per node, n / p = ", ratios[j], "\n")
  for (nrel in 1:N_realizations) {
    E(b1)$weight <- runif(gsize(b1), min = 0, max = 1)
    Ltrue <- as.matrix(laplacian_matrix(b1))
    Atrue <- diag(diag(Ltrue)) - Ltrue
    Y <- MASS::mvrnorm(T, mu = rep(0, n1 + n2), Sigma = MASS::ginv(Ltrue))
    S <- cov(Y)
    for (i in gamma_size) {
      print(gamma_set[i])
      graph <- learn_bipartite_graph(S, w0 = "qp", z = n1-n2, nu = gamma_set[i], maxiter = 100000)
      gamma_rel_err[i, j] <- gamma_rel_err[i, j] + relativeError(Ltrue, graph$Laplacian)
      gamma_fscore[i, j] <- gamma_fscore[i, j] + metrics(Ltrue, graph$Laplacian, 5e-2)[1]
    }
  }
  gamma_rel_err[,j] <- gamma_rel_err[,j] / N_realizations
  gamma_fscore[,j] <- gamma_fscore[,j] / N_realizations
}

saveRDS(gamma_rel_err, file = "gamma_rel_err.rds")
saveRDS(gamma_fscore, file = "gamma_fscore.rds")

#colors <- c("#0B032D", "#843B62", "#F67E7D", "#FFB997")
#gr = .5 * (1 + sqrt(5))
#setEPS()
#postscript("relative_error_gamma.ps", family = "ComputerModern", height = 5.5, width = gr * 4.5, pointsize = 14)
#plot(gamma_size, gamma_rel_err[,1], type = "b", pch=15, cex=.75, ylim=c(0, max(gamma_rel_err)),
#     xlab = TeX("$\\gamma$"), ylab = "Average Relative Error", col = colors[1], xaxt = "n")
#grid()
#lines(gamma_size, gamma_rel_err[,2], type = "b", pch=16, cex=.75, col = colors[2], xaxt = "n")
#lines(gamma_size, gamma_rel_err[,3], type = "b", pch=17, cex=.75, col = colors[3], xaxt = "n")
#lines(gamma_size, gamma_rel_err[,4], type = "b", pch=18, cex=.85, col = colors[4], xaxt = "n")
#axis(side = 1, at = gamma_size, labels = gamma_set)
#legend("topright", legend = c(TeX("$n / p = .5$"), TeX("$n / p = 10$"), TeX("$n / p = 100$"), TeX("$n / p = 1000$")),
#       col=colors, pch=c(15, 16, 17, 18), lty=c(1, 1, 1, 1), bty="n")
#dev.off()
#embed_fonts("relative_error_gamma.ps", outfile="relative_error_gamma.ps")
#setEPS()
#postscript("fscore_gamma.ps", family = "ComputerModern", height = 5.5, width = gr * 4.5, pointsize = 14)
#plot(gamma_size, gamma_fscore[,1], ylim=c(min(gamma_fscore), 1.), xlab = TeX("$\\gamma$"),
#     ylab = "Average F-score", type = "b", pch=15, cex=.75, col = colors[1], xaxt = "n")
#grid()
#lines(gamma_size, gamma_fscore[,2], type = "b", pch=16, cex=.75, col = colors[2], xaxt = "n")
#lines(gamma_size, gamma_fscore[,3], type = "b", pch=17, cex=.75, col = colors[3], xaxt = "n")
#lines(gamma_size, gamma_fscore[,4], type = "b", pch=18, cex=.85, col = colors[4], xaxt = "n")
#axis(side = 1, at = gamma_size, labels = gamma_set)
#legend("bottomright", legend = c(TeX("$n / p = .5$"), TeX("$n / p = 10$"), TeX("$n / p = 100$"), TeX("$n / p = 1000$")),
#       col=colors, pch=c(15, 16, 17, 18), lty=c(1, 1, 1, 1), bty="n")
#dev.off()
#embed_fonts("fscore_gamma.ps", outfile="fscore_gamma.ps")
