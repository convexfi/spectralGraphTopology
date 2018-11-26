library(igraph)
library(viridis)
library(corrplot)
library(spectralGraphTopology)
library(glasso)
library(latex2exp)

N_realizations <- 100
N <- 64
T <- 30 * N
K <- 4
P <- diag(1, K)
# probability of node connection
p <- 1
kappa_seq <- c(.6, .7, .8)
#kappa_seq <- c(.1, .2, .3, .4, .5, .6, .7, .8)
# K-component graph
mgraph <- sample_sbm(N, pref.matrix = P, block.sizes = c(rep(N / K, K)))
erdos_renyi <- erdos.renyi.game(N, p)
rel_err <- array(0, length(kappa_seq))
fscore <- array(0, length(kappa_seq))
for (k in c(1:length(kappa_seq))) {
  a <- kappa_seq[k]
  for (n in 1:N_realizations) {
    cat("\nRunning simulation for kappa =", a, "\n")
    E(mgraph)$weight <- runif(gsize(mgraph), min = 0, max = 1)
    Ltrue <- as.matrix(laplacian_matrix(mgraph))
    Wtrue <- diag(diag(Ltrue)) - Ltrue
    # Erdo-Renyi as noise model
    E(erdos_renyi)$weight <- runif(gsize(erdos_renyi), min = 0, max = a)
    Lerdo <- as.matrix(laplacian_matrix(erdos_renyi))
    # Noisy Laplacian
    Lnoisy <- Ltrue + Lerdo
    Y <- MASS::mvrnorm(T, mu = rep(0, N), Sigma = MASS::ginv(Lnoisy))
    graph <- learnGraphTopology(crossprod(Y) / T, w0 = "naive", K = 4, beta = 4)
    re <- relativeError(Ltrue, graph$Lw)
    fs <- Fscore(Ltrue, graph$Lw, 1e-3)
    rel_err[k] <- rel_err[k] + re
    fscore[k] <- fscore[k] + fs
    print(re)
    print(fs)
    print(graph$convergence)
  }
  rel_err[k] <- rel_err[k] / N_realizations
  fscore[k] <- fscore[k] / N_realizations
  cat("\n** spectralGraphTopology results **\n")
  cat("Avg Relative error: ")
  cat(rel_err[k], "\n")
  cat("Avg Fscore: ")
  cat(fscore[k], "\n")
}
gr = .5 * (1 + sqrt(5))
setEPS()
postscript("relative_error.ps", family = "Times", height = 5, width = gr * 3.5)
plot(c(1:length(kappa_seq)), rel_err, type = "b", pch=19, cex=.75,
     ylim=c(min(rel_err) - 1e-2, max(rel_err) + .1e-2),
     xlab = TeX("$\\kappa$"), ylab = "Average Relative Error", col = "black", xaxt = "n")
axis(side = 1, at = c(1:length(kappa_seq)), labels = kappa_seq)
dev.off()
setEPS()
postscript("fscore.ps", family = "Times", height = 5, width = gr * 3.5)
plot(c(1:length(kappa_seq)), fscore, type = "b", pch=19, cex=.75,
     ylim=c(min(fscore) - 1e-2, 1.),
     xlab = TeX("$\\kappa$"), ylab = "Average F-score", col = "black", xaxt = "n")
axis(side = 1, at = c(1:length(kappa_seq)), labels = kappa_seq)
dev.off()
