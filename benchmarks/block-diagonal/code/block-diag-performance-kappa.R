library(igraph)
library(viridis)
library(corrplot)
library(spectralGraphTopology)
library(latex2exp)
library(extrafont)
set.seed(42)

N_realizations <- 10
N <- 20
T <- 30 * N
K <- 4
P <- diag(1, K)
# probability of node connection
kappa_seq <- c(.1, .2, .3, .4, .5, .6)
# K-component graph
mgraph <- sample_sbm(N, pref.matrix = P, block.sizes = c(rep(N / K, K)))
rel_err <- array(0, length(kappa_seq))
fscore <- array(0, length(kappa_seq))
for (k in c(1:length(kappa_seq))) {
  a <- kappa_seq[k]
  cat("\nRunning simulation for kappa =", a, "\n")
  for (n in 1:N_realizations) {
    E(mgraph)$weight <- runif(gsize(mgraph), min = 0, max = 1)
    Ltrue <- as.matrix(laplacian_matrix(mgraph))
    # Erdo-Renyi as noise model
    erdos_renyi <- erdos.renyi.game(N, p = .35)
    E(erdos_renyi)$weight <- runif(gsize(erdos_renyi), min = 0, max = a)
    Lerdo <- as.matrix(laplacian_matrix(erdos_renyi))
    # Noisy Laplacian
    Lnoisy <- Ltrue + Lerdo
    Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Lnoisy))
    S <- cov(Y)
    graph <- learn_laplacian_matrix(S, w0 = "qp", k = 4, beta = 400, fix_beta = TRUE,
                                    maxiter = 100000, abstol = 0)
    print(graph$convergence)
    fs <- metrics(Ltrue, graph$Laplacian, 1e-3)[1]
    print(fs)
    re <- relativeError(Ltrue, graph$Laplacian)
    rel_err[k] <- rel_err[k] + re
    fscore[k] <- fscore[k] + fs
  }
  rel_err[k] <- rel_err[k] / N_realizations
  fscore[k] <- fscore[k] / N_realizations
  cat("\n** spectralGraphTopology results **\n")
  cat("Avg Relative error: ")
  cat(rel_err[k], "\n")
  cat("Avg Fscore: ")
  cat(fscore[k], "\n")
}
gr <- .5 * (1 + sqrt(5))
legend <- c("Average Relative Error", "Average F-score")
colors <- c("#ff5252", "black")
setEPS()
postscript("../latex/figures/performance_kappa.ps", family = "ComputerModern", height = 5, width = gr * 4)
par(mar = c(5, 5, 3, 5))
plot(c(1:length(kappa_seq)), rel_err, type = "b", pch=19, cex=.75, ylim=c(min(rel_err), max(rel_err)),
     xlab = TeX("$\\kappa$"), ylab = "Average Relative Error", col = colors[1], xaxt = "n", lty = 1, lwd = 2)
axis(side = 1, at = c(1:length(kappa_seq)), labels = kappa_seq)
par(new = TRUE)
plot(c(1:length(kappa_seq)), fscore, type = "b", pch=19, cex=.75, ylim=c(min(fscore), 1.),
     xlab = "", ylab = "", col = colors[2], xaxt = "n", lwd = 2, xaxt = "n", yaxt = "n", lty = 3)
grid()
axis(side = 4)
mtext("Average F-score", side = 4, line = 3)
legend("top", legend = legend, col = colors, lty = c(1, 3), bty = "n", lwd = c(2, 2))
dev.off()
embed_fonts("../latex/figures/performance_kappa.ps", outfile="../latex/figures/performance_kappa.ps")
